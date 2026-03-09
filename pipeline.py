#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import subprocess
import sys

try:
	from bids import BIDSLayout
except Exception:
	BIDSLayout = None

try:
	from utils import app_name, version, release_date
except Exception:
	app_name = 'GGR-recon'
	version = 'unknown'
	release_date = ''

ACQ_ORDER = ['sag', 'cor', 'ax']
FILTER_KEY_ALIASES = {
	'sub': 'subject',
	'ses': 'session',
	'acq': 'acquisition',
	'rec': 'reconstruction',
}
GROUP_EXCLUDED_ENTITIES = {'acquisition', 'suffix', 'extension', 'datatype'}

def print_help():
	print('usage: pipeline.py [PREPROCESS_ARGS ...] [-- RECON_ARGS ...]')
	print('')
	print('Runs preprocess.py first, then recon.py.')
	print('Arguments before "--" are passed to preprocess.py.')
	print('Arguments after "--" are passed to recon.py.')
	print('Pipeline-owned options:')
	print('  --participant-label LABEL [LABEL ...]  map to --bids-filter subject=...')
	print('  --session-label LABEL [LABEL ...]      map to --bids-filter session=...')
	print('  --bids-filter-file FILE                nested BIDS filter JSON with a top-level "t2w" block')
	print('')
	print('Examples:')
	print('  pipeline.py --path /data --temp_path /temp --out_path /bids')
	print('  pipeline.py --path /data --temp_path /temp --out_path /bids \\')
	print('    --bids-filter subject=2983 --bids-filter rec=filtered -- --ggr -w 0.03')
	print('  pipeline.py --path /data --participant-label 2983 --session-label 1a -- --ggr -w 0.03')
	print('')
	print('Notes:')
	print('  - If "--" is omitted, no extra args are passed to recon.py (defaults are used).')
	print('  - If no explicit -f/--filenames is provided, pipeline runs all complete BIDS groups matching filters.')
	print('  - In --bids-filter-file conflicts, file values override --bids-filter values.')
	print('  - participant/session labels override subject/session values from --bids-filter.')
	print('  - In -f/--filenames mode, participant/session labels and --bids-filter-file are ignored with a warning.')
	print('  - All original preprocess.py and recon.py arguments are supported via passthrough.')


def split_passthrough_args(argv):
	if '--' in argv:
		sep = argv.index('--')
		return argv[:sep], argv[sep + 1:]
	return argv, []

def flatten_label_values(values):
	flat = []
	for chunk in values or []:
		items = chunk if isinstance(chunk, (list, tuple)) else [chunk]
		for item in items:
			value = str(item).strip()
			if value != '' and value not in flat:
				flat.append(value)
	return flat

def parse_pipeline_options(preprocess_args):
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument('--participant-label', action='append', nargs='+', default=[])
	parser.add_argument('--session-label', action='append', nargs='+', default=[])
	parser.add_argument('--bids-filter-file')
	parsed, remaining = parser.parse_known_args(preprocess_args)
	participants = flatten_label_values(parsed.participant_label)
	sessions = flatten_label_values(parsed.session_label)
	return remaining, participants, sessions, parsed.bids_filter_file

def normalize_filter_key(key):
	key = str(key).strip()
	if key == '':
		return None
	return FILTER_KEY_ALIASES.get(key, key)

def normalize_filter_value(value):
	if isinstance(value, list):
		items = []
		for item in value:
			if isinstance(item, (dict, list)) or item is None:
				raise ValueError('invalid filter value type in --bids-filter-file: %s' % type(item).__name__)
			text = str(item).strip()
			if text != '':
				items.append(text)
		if len(items) == 0:
			raise ValueError('empty list value in --bids-filter-file is not allowed')
		return items
	if isinstance(value, dict) or value is None:
		raise ValueError('invalid filter value type in --bids-filter-file: %s' % type(value).__name__)
	text = str(value).strip()
	if text == '':
		raise ValueError('empty filter value in --bids-filter-file is not allowed')
	return text

def load_bids_filters_from_file(path):
	try:
		with open(path, 'r') as f:
			data = json.load(f)
	except OSError as exc:
		raise ValueError('could not read --bids-filter-file "%s": %s' % (path, str(exc)))
	except ValueError as exc:
		raise ValueError('invalid JSON in --bids-filter-file "%s": %s' % (path, str(exc)))

	if not isinstance(data, dict):
		raise ValueError('--bids-filter-file must contain a JSON object')
	if 't2w' not in data:
		raise ValueError('--bids-filter-file must contain a top-level "t2w" object')
	block = data.get('t2w')
	if not isinstance(block, dict):
		raise ValueError('the "t2w" value in --bids-filter-file must be a JSON object')

	parsed = {}
	for raw_key, raw_value in block.items():
		key = normalize_filter_key(raw_key)
		if key is None:
			raise ValueError('invalid empty key in --bids-filter-file')
		parsed[key] = normalize_filter_value(raw_value)
	return parsed

def parse_preprocess_path(args):
	path = '/opt/GGR-recon/data/'
	ii = 0
	while ii < len(args):
		token = args[ii]
		if token in ('-p', '--path') and ii + 1 < len(args):
			path = args[ii + 1]
			ii += 2
			continue
		if token.startswith('--path='):
			path = token.split('=', 1)[1]
		ii += 1
	return path

def has_filenames_arg(args):
	return '-f' in args or '--filenames' in args

def has_option(args, names):
	for ii, token in enumerate(args):
		if token in names:
			return True
		for name in names:
			if token.startswith(name + '='):
				return True
	return False

def get_last_option_value(args, names):
	value = None
	ii = 0
	while ii < len(args):
		token = args[ii]
		if token in names and ii + 1 < len(args):
			value = args[ii + 1]
			ii += 2
			continue
		for name in names:
			if token.startswith(name + '='):
				value = token.split('=', 1)[1]
				break
		ii += 1
	return value

def extract_bids_filters(args):
	raw_filters = []
	ii = 0
	while ii < len(args):
		token = args[ii]
		if token == '--bids-filter' and ii + 1 < len(args):
			raw_filters.append(args[ii + 1])
			ii += 2
			continue
		if token.startswith('--bids-filter='):
			raw_filters.append(token.split('=', 1)[1])
		ii += 1
	return raw_filters

def parse_filter_key_value(raw):
	if '=' not in raw:
		return None, None
	key, value = raw.split('=', 1)
	key = FILTER_KEY_ALIASES.get(key.strip(), key.strip())
	value = value.strip()
	if key == '' or value == '':
		return None, None
	if ',' in value:
		value = [v.strip() for v in value.split(',') if v.strip() != '']
	return key, value

def remove_bids_filter_keys(args, keys_to_remove):
	out = []
	ii = 0
	while ii < len(args):
		token = args[ii]
		if token == '--bids-filter':
			if ii + 1 < len(args):
				raw = args[ii + 1]
				key, _ = parse_filter_key_value(raw)
				if key in keys_to_remove:
					ii += 2
					continue
				out += [token, raw]
				ii += 2
				continue
			out.append(token)
			ii += 1
			continue
		if token.startswith('--bids-filter='):
			raw = token.split('=', 1)[1]
			key, _ = parse_filter_key_value(raw)
			if key in keys_to_remove:
				ii += 1
				continue
		out.append(token)
		ii += 1
	return out

def apply_label_filters(preprocess_args, participant_labels, session_labels):
	args = list(preprocess_args)
	remove_keys = set()
	if len(participant_labels) > 0:
		remove_keys.add('subject')
	if len(session_labels) > 0:
		remove_keys.add('session')
	if len(remove_keys) > 0:
		args = remove_bids_filter_keys(args, remove_keys)
	if len(participant_labels) > 0:
		args += ['--bids-filter', 'subject=%s' % ','.join(participant_labels)]
	if len(session_labels) > 0:
		args += ['--bids-filter', 'session=%s' % ','.join(session_labels)]
	return args

def apply_file_filters(preprocess_args, file_filters):
	args = list(preprocess_args)
	if len(file_filters) == 0:
		return args
	args = remove_bids_filter_keys(args, set(file_filters.keys()))
	for key, value in file_filters.items():
		if isinstance(value, list):
			args += ['--bids-filter', '%s=%s' % (key, ','.join(value))]
		else:
			args += ['--bids-filter', '%s=%s' % (key, value)]
	return args

def group_key_from_entities(entities):
	items = []
	for key, value in entities.items():
		if value is None or key in GROUP_EXCLUDED_ENTITIES:
			continue
		items.append((key, str(value)))
	return tuple(sorted(items))

def better_path(path_a, path_b):
	if path_a is None:
		return path_b
	depth_a = path_a.count(os.sep)
	depth_b = path_b.count(os.sep)
	if depth_b < depth_a:
		return path_b
	if depth_a < depth_b:
		return path_a
	return min(path_a, path_b)

def format_group_key(group_key):
	order = {'subject': 0, 'session': 1, 'reconstruction': 2}
	pairs = sorted(group_key, key=lambda kv: (order.get(kv[0], 99), kv[0], kv[1]))
	return '_'.join('%s-%s' % (key, value) for key, value in pairs)

def discover_group_filter_sets(preprocess_args):
	if BIDSLayout is None:
		return None

	root = parse_preprocess_path(preprocess_args)
	try:
		layout = BIDSLayout(root, validate=False)
	except Exception:
		return []

	query = {
		'suffix': 'T2w',
		'acquisition': ACQ_ORDER,
		'extension': ['.nii', '.nii.gz'],
		'datatype': 'anat',
		'scope': 'raw',
	}
	for raw_filter in extract_bids_filters(preprocess_args):
		key, value = parse_filter_key_value(raw_filter)
		if key is not None:
			query[key] = value

	try:
		bids_files = layout.get(return_type='object', **query)
	except Exception:
		return []
	groups = {}
	for bids_file in bids_files:
		entities = bids_file.get_entities()
		acq = str(entities.get('acquisition', ''))
		if acq not in ACQ_ORDER:
			continue
		if entities.get('subject') is None:
			continue

		group_key = group_key_from_entities(entities)
		if group_key not in groups:
			groups[group_key] = {'acq_map': {}}
		current = groups[group_key]['acq_map'].get(acq)
		groups[group_key]['acq_map'][acq] = better_path(current, bids_file.path)

	complete = []
	for group_key, group in groups.items():
		if all(acq in group['acq_map'] for acq in ACQ_ORDER):
			filter_args = []
			for key, value in group_key:
				filter_args += ['--bids-filter', '%s=%s' % (key, value)]
			complete.append((group_key, filter_args))

	complete.sort(key=lambda item: format_group_key(item[0]))
	return complete

def run_single(preprocess_args, recon_args):
	final_recon_args = list(recon_args)
	preprocess_temp_names = ['-t', '--temp_path', '--working_path', '-w']
	recon_temp_names = ['-t', '--temp_path', '--working_path']
	out_names = ['-o', '--out_path']

	if not has_option(final_recon_args, recon_temp_names):
		temp_value = get_last_option_value(preprocess_args, preprocess_temp_names)
		if temp_value is not None:
			final_recon_args = ['--temp_path', temp_value] + final_recon_args

	if not has_option(final_recon_args, out_names):
		out_value = get_last_option_value(preprocess_args, out_names)
		if out_value is not None:
			final_recon_args = ['--out_path', out_value] + final_recon_args

	rc = run_script('preprocess.py', preprocess_args)
	if rc != 0:
		print('[pipeline] preprocess.py failed with exit code %d' % rc)
		return rc

	rc = run_script('recon.py', final_recon_args)
	if rc != 0:
		print('[pipeline] recon.py failed with exit code %d' % rc)
		return rc
	return 0

def run_script(script_name, script_args):
	script_path = os.path.join(os.path.dirname(__file__), script_name)
	cmd = [sys.executable, script_path] + script_args
	print('[pipeline] running:', ' '.join(shlex.quote(token) for token in cmd))
	result = subprocess.run(cmd)
	return result.returncode


def main():
	argv = sys.argv[1:]

	if '-h' in argv or '--help' in argv:
		print_help()
		return 0
	if '-V' in argv or '--version' in argv:
		print('%s version : v %s %s' % (app_name, version, release_date))
		return 0

	preprocess_args, recon_args = split_passthrough_args(argv)
	preprocess_args, participant_labels, session_labels, bids_filter_file = parse_pipeline_options(preprocess_args)
	if has_filenames_arg(preprocess_args):
		if len(participant_labels) > 0 or len(session_labels) > 0:
			print('[pipeline] warning: --participant-label/--session-label are ignored when -f/--filenames is used.')
		if bids_filter_file is not None:
			print('[pipeline] warning: --bids-filter-file is ignored when -f/--filenames is used.')
	else:
		if bids_filter_file is not None:
			try:
				file_filters = load_bids_filters_from_file(bids_filter_file)
			except ValueError as exc:
				print('Error:', str(exc))
				return 1
			preprocess_args = apply_file_filters(preprocess_args, file_filters)
		preprocess_args = apply_label_filters(preprocess_args, participant_labels, session_labels)
	# Expand into all matching groups unless explicit filenames are provided.
	# This includes cases with filters (e.g., subject/session without rec).
	should_expand_groups = not has_filenames_arg(preprocess_args)

	if should_expand_groups:
		discovered = discover_group_filter_sets(preprocess_args)
		if discovered is None:
			print('[pipeline] pybids is unavailable; running a single preprocess/recon pair.')
			return run_single(preprocess_args, recon_args)
		if len(discovered) == 0:
			print('[pipeline] no complete BIDS groups found with no filters; running single pass.')
			return run_single(preprocess_args, recon_args)

		print('[pipeline] discovered %d complete BIDS groups.' % len(discovered))
		for idx, (group_key, group_filter_args) in enumerate(discovered, start=1):
			print('[pipeline] group %d/%d: %s' % (
					idx, len(discovered), format_group_key(group_key)))
			rc = run_single(preprocess_args + group_filter_args, recon_args)
			if rc != 0:
				return rc
		return 0

	return run_single(preprocess_args, recon_args)


if __name__ == '__main__':
	sys.exit(main())
