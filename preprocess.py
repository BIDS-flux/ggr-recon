#!/usr/bin/env python3

import numpy as np
import SimpleITK as sitk
from scipy.io import loadmat, savemat
from scipy import signal
import json
import os
import sys
import time
import argparse
import pathlib

try:
	from bids import BIDSLayout
	from bids.layout import parse_file_entities
except ImportError:
	BIDSLayout = None
	parse_file_entities = None

from utils import *

from rich.console import Console
from rich import box
from rich.table import Table
from rich.panel import Panel
from rich.progress import track

ACQ_ORDER = ['sag', 'cor', 'ax']
FILTER_KEY_ALIASES = {
	'sub': 'subject',
	'ses': 'session',
	'acq': 'acquisition',
	'rec': 'reconstruction',
}
ENTITY_NAME_TO_LABEL = {
	'subject': 'sub',
	'session': 'ses',
	'acquisition': 'acq',
	'reconstruction': 'rec',
	'run': 'run',
	'echo': 'echo',
	'part': 'part',
	'desc': 'desc',
	'space': 'space',
	'ceagent': 'ce',
	'direction': 'dir',
	'chunk': 'chunk',
	'task': 'task',
	'inversion': 'inv',
}
KNOWN_DATATYPES = {
	'anat', 'func', 'dwi', 'fmap', 'perf', 'pet', 'meg', 'eeg', 'ieeg', 'micr'
}
GROUP_EXCLUDED_ENTITIES = {'acquisition', 'suffix', 'extension'}

def ensure_dir(path):
	if path.endswith('/'):
		return path
	return path + '/'

def normalize_bids_entity(name, value):
	if value is None:
		return None
	value = str(value)
	prefix = name + '-'
	if value.startswith(prefix):
		return value
	return prefix + value

def entity_to_token(key, value):
	label = ENTITY_NAME_TO_LABEL.get(key, key)
	return '%s-%s' % (label, value)

def infer_datatype_from_path(filename):
	parent = pathlib.Path(filename).parent.name
	if parent in KNOWN_DATATYPES:
		return parent
	for part in reversed(pathlib.Path(filename).parts):
		if part in KNOWN_DATATYPES:
			return part
	return 'anat'

def parse_bids_filters(raw_filters):
	parsed = {}
	for raw_filter in raw_filters:
		if '=' not in raw_filter:
			raise ValueError('invalid --bids-filter "%s" (expected KEY=VALUE)' % raw_filter)
		key, value = raw_filter.split('=', 1)
		key = key.strip()
		value = value.strip()
		if key == '' or value == '':
			raise ValueError('invalid --bids-filter "%s" (expected KEY=VALUE)' % raw_filter)
		key = FILTER_KEY_ALIASES.get(key, key)
		if ',' in value:
			parsed[key] = [v.strip() for v in value.split(',') if v.strip() != '']
		else:
			parsed[key] = value
	return parsed

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

def group_key_from_entities(entities):
	items = []
	for key, value in entities.items():
		if value is None or key in GROUP_EXCLUDED_ENTITIES:
			continue
		items.append((key, str(value)))
	return tuple(sorted(items))

def format_group_label(group_entities):
	parts = []
	sub = normalize_bids_entity('sub', group_entities.get('subject'))
	ses = normalize_bids_entity('ses', group_entities.get('session'))
	if sub is not None:
		parts.append(sub)
	if ses is not None:
		parts.append(ses)

	for key in sorted(group_entities.keys()):
		if key in {'subject', 'session', 'datatype'}:
			continue
		if key in GROUP_EXCLUDED_ENTITIES:
			continue
		value = group_entities.get(key)
		if value is None:
			continue
		parts.append(entity_to_token(key, value))

	return '_'.join(parts)

def collect_candidate_groups(file_records):
	groups = {}
	for record in file_records:
		entities = record['entities']
		acq = str(entities.get('acquisition', ''))
		if acq not in ACQ_ORDER:
			continue

		group_key = group_key_from_entities(entities)
		if group_key not in groups:
			groups[group_key] = {'entities': dict(entities), 'acq_map': {}}

		if acq in groups[group_key]['acq_map']:
			prev = groups[group_key]['acq_map'][acq]
			raise ValueError('duplicate "%s" acquisition found for %s: %s and %s' % (
					acq, format_group_label(groups[group_key]['entities']), prev, record['path']))
		groups[group_key]['acq_map'][acq] = record['path']

	return groups

def group_priority(group):
	score = 0
	if group['entities'].get('datatype') == 'anat':
		score += 4
	paths = sorted(group['acq_map'].values())
	# Prefer shallower paths when equivalent groups are duplicated.
	score -= int(sum(path.count(os.sep) for path in paths) / max(1, len(paths)))
	return score

def choose_complete_group(groups):
	complete_groups = []
	for group in groups.values():
		if all(acq in group['acq_map'] for acq in ACQ_ORDER):
			complete_groups.append(group)

	# Collapse duplicate logical groups that differ only by index context/path.
	dedup_complete = {}
	for group in complete_groups:
		label = format_group_label(group['entities'])
		existing = dedup_complete.get(label)
		if existing is None or group_priority(group) > group_priority(existing):
			dedup_complete[label] = group
	complete_groups = list(dedup_complete.values())

	if len(complete_groups) == 0:
		raise ValueError('no complete acq-{sag,cor,ax} set found for suffix T2w')
	if len(complete_groups) > 1:
		labels = sorted([format_group_label(group['entities']) for group in complete_groups])
		raise ValueError('multiple complete BIDS groups found: %s. Use --bids-filter (e.g., rec=filtered or rec=orig) to select one.' % ', '.join(labels))

	return complete_groups[0]

def build_bids_output_name(reference_file, source_reconstruction=None):
	name = os.path.basename(reference_file)
	if name.endswith('.nii.gz'):
		ext = '.nii.gz'
		stem = name[:-7]
	elif name.endswith('.nii'):
		ext = '.nii'
		stem = name[:-4]
	else:
		return None
	tokens = stem.split('_')
	if len(tokens) < 2 or tokens[-1] != 'T2w':
		return None

	rec_value = None
	if source_reconstruction is not None:
		rec_value = str(source_reconstruction)

	out_tokens = []
	for token in tokens[:-1]:
		if token.startswith('acq-') or token.startswith('rec-'):
			if token.startswith('rec-') and rec_value is None:
				rec_value = token.split('-', 1)[1]
			continue
		out_tokens.append(token)
	if rec_value is not None and rec_value != '':
		# Preserve source rec-* information without conflicting with rec-superesolution.
		out_tokens.append('acq-' + rec_value)
	out_tokens.append('rec-superesolution')
	out_tokens.append('T2w')
	return '_'.join(out_tokens) + ext

def relativize_paths(paths, root_path):
	if root_path is None:
		return [os.path.abspath(path) for path in paths]

	output = []
	for path in paths:
		abspath = os.path.abspath(path)
		try:
			relpath = os.path.relpath(abspath, os.path.abspath(root_path))
			if not relpath.startswith('..'):
				output.append(relpath.replace(os.sep, '/'))
			else:
				output.append(abspath)
		except ValueError:
			output.append(abspath)
	return output

def build_bids_info(group, flist, root_path=None):
	entities = group['entities']
	sub = normalize_bids_entity('sub', entities.get('subject'))
	ses = normalize_bids_entity('ses', entities.get('session'))
	if sub is None:
		return None

	output_rel_dir = os.path.join(sub, 'anat')
	if ses is not None:
		output_rel_dir = os.path.join(sub, ses, 'anat')

	output_name = build_bids_output_name(
			group['acq_map'].get('sag', flist[0]),
			source_reconstruction=entities.get('reconstruction'))
	if output_name is None:
		if ses is None:
			output_name = '%s_rec-superesolution_T2w.nii.gz' % sub
		else:
			output_name = '%s_%s_rec-superesolution_T2w.nii.gz' % (sub, ses)

	relevant_entities = {}
	for key, value in entities.items():
		if key in GROUP_EXCLUDED_ENTITIES or value is None:
			continue
		if key == 'datatype':
			continue
		relevant_entities[key] = str(value)

	return {
		'output_name': output_name,
		'output_rel_dir': output_rel_dir.replace(os.sep, '/'),
		'subject': sub,
		'session': ses,
		'datatype': 'anat',
		'input_acquisitions': ACQ_ORDER,
		'source_images': relativize_paths(flist, root_path),
		'source_entities': relevant_entities,
	}

def discover_bids_inputs(search_path, extra_filters):
	if BIDSLayout is None:
		raise ImportError('pybids is required for automatic BIDS discovery. Install "pybids".')

	layout = BIDSLayout(search_path, validate=False)
	filters = {
		'suffix': 'T2w',
		'extension': ['.nii', '.nii.gz'],
		'scope': 'raw',
	}
	if 'acquisition' not in extra_filters:
		filters['acquisition'] = ACQ_ORDER
	if 'datatype' not in extra_filters:
		filters['datatype'] = 'anat'
	filters.update(extra_filters)

	bids_files = layout.get(return_type='object', **filters)
	file_records = []
	for bids_file in bids_files:
		entities = bids_file.get_entities()
		if entities.get('subject') is None:
			continue
		file_records.append({'path': bids_file.path, 'entities': entities})

	groups = collect_candidate_groups(file_records)
	selected_group = choose_complete_group(groups)
	flist = [selected_group['acq_map'][acq] for acq in ACQ_ORDER]
	bids_info = build_bids_info(selected_group, flist, root_path=search_path)
	return flist, bids_info

def select_bids_inputs_from_filenames(filenames):
	if parse_file_entities is None:
		raise ImportError('pybids is required for BIDS entity parsing. Install "pybids".')

	file_records = []
	for filename in filenames:
		entities = parse_file_entities(filename)
		if entities.get('suffix') != 'T2w':
			continue
		if str(entities.get('acquisition', '')) not in ACQ_ORDER:
			continue
		if entities.get('subject') is None:
			continue
		entities = dict(entities)
		entities['datatype'] = infer_datatype_from_path(filename)
		file_records.append({'path': filename, 'entities': entities})

	if len(file_records) == 0:
		return filenames, None
	if len(file_records) != len(filenames):
		raise ValueError('explicit filenames must all be BIDS-compatible when using acq-{sag,cor,ax}_T2w input selection')

	groups = collect_candidate_groups(file_records)
	selected_group = choose_complete_group(groups)
	flist = [selected_group['acq_map'][acq] for acq in ACQ_ORDER]

	common_root = os.path.commonpath([os.path.abspath(path) for path in flist])
	bids_info = build_bids_info(selected_group, flist, root_path=common_root)
	return flist, bids_info

parser = argparse.ArgumentParser()
parser.add_argument('-V', '--version', action='version',
		version='%s version : v %s %s' % (app_name, version, release_date),
		help='show version')

parser.add_argument('-f', '--filenames', nargs='+',
        help='filenames of input the low-res images; (full path required)\
                e.g., -f a.nii.gz b.nii.gz c.nii.gz')
parser.add_argument('-s', '--size', nargs='+', type=int,
		help='size of the high-res reconstruction, optional; \
				even positive integers required if set; \
                e.g., -s 312 384 330')
parser.add_argument('-r', '--resample', action='store_true',
		help='resample the first low-res image in the high-res lattice \
				and then exit. Usually used for determining a user \
				defined size of the high-res reconstruction')
parser.add_argument('-p', '--path', default='/opt/GGR-recon/data/')
parser.add_argument('-t', '--temp_path', '-w', '--working_path',
		default='/opt/GGR-recon/temp/',
		help='path for intermediate files; default is /opt/GGR-recon/temp/')
parser.add_argument('-o', '--out_path', default='/opt/GGR-recon/recons/')
parser.add_argument('--bids-filter', action='append', default=[],
		help='additional pybids filters for automatic discovery, as KEY=VALUE (repeatable)')
parser.add_argument('--bids-filter-file',
		help='path to nested JSON filters (expects top-level "t2w" object)')
args = parser.parse_args()
flist = args.filenames
sz = args.size
resample_only = args.resample

path = ensure_dir(args.path)
working_path = ensure_dir(args.temp_path)
out_path = ensure_dir(args.out_path)

bids_filters = {}
if flist is not None and len(flist) > 0 and args.bids_filter_file is not None:
	print('Warning: --bids-filter-file is ignored when -f/--filenames is used.')
else:
	try:
		bids_filters = parse_bids_filters(args.bids_filter)
		if args.bids_filter_file is not None:
			file_filters = load_bids_filters_from_file(args.bids_filter_file)
			bids_filters.update(file_filters)
	except ValueError as exc:
		print('Error:', str(exc))
		sys.exit(1)

bids_info = None
if flist is None or len(flist) == 0:
	try:
		flist, bids_info = discover_bids_inputs(path, bids_filters)
	except (ImportError, ValueError) as exc:
		print('Error:', str(exc))
		sys.exit(1)
	if flist is None:
		print('No BIDS image data found at:', path)
		print('Expected at least one complete acq-{sag,cor,ax}_T2w set per subject[/session]')
		print('Use --bids-filter to disambiguate if multiple sets exist')
		sys.exit(1)
else:
	try:
		flist, bids_info = select_bids_inputs_from_filenames(flist)
	except (ImportError, ValueError) as exc:
		print('Error:', str(exc))
		sys.exit(1)

n_imgs = len(flist)
if n_imgs == 0:
	print('No image data found!')
	sys.exit(1)

if sz != None and (len(sz) != 3 or np.any(np.array(sz) <= 0)):
	print('SIZE =', sz)
	print('Error: SIZE should comprise 3 positive integers')
	sys.exit(1)

print('path : ' + str(path))
print('temp_path : ' + str(working_path))
print('out_path : ' + str(out_path))

os.makedirs(out_path, exist_ok=True)
os.makedirs(working_path, exist_ok=True)

bids_output_file = os.path.join(working_path, 'bids_output_name.json')
if bids_info is not None:
	with open(bids_output_file, 'w') as f:
		json.dump(bids_info, f, indent=2)
else:
	if os.path.isfile(bids_output_file):
		os.remove(bids_output_file)

img_path = []
img_fn = []
img_ext = []


for filename in flist:
    fpname =  pathlib.PurePosixPath(filename)
    base, first_dot, rest = fpname.name.partition('.')
    #filename = filename.with_name(base)
    p = str(fpname.parent)
    if not p.endswith('/'):
        p += '/'
    img_path.append( p )
    img_fn.append( str(base) )
    img_ext.append( '.'+str(rest) )

console = Console()
print_header(console)

# step 0: make the orientations the same for all LR images
for ii in range(0, n_imgs):
  inputVolume = img_path[ii] + img_fn[ii] + img_ext[ii]
  outputVolume = working_path + img_fn[ii] + img_ext[ii]
  print(str(inputVolume))
  print(str(outputVolume))
  reader = sitk.ImageFileReader()
  reader.SetFileName( inputVolume )
  inputImage = reader.Execute();
  # Now we clone the input image.
  reorientedImage = sitk.DICOMOrient( inputImage, 'LPS' )
  writer = sitk.ImageFileWriter()
  writer.SetFileName( outputVolume )
  writer.Execute( reorientedImage )

#print('completed step 0')
#print('\t- make the orientations the same for all LR images')


# step 1: resample the images
img0 = imread(working_path + img_fn[0] + img_ext[0])
if sz == None:
	img0x = resample_iso_img(img0)
	sz = img0x.GetSize()
else:
	img0x = resample_iso_img_with_size(img0, sz)

sz = img0x.GetSize() # update the variable of image size

# =========== Print summary of the execution =============
mode = 'Preprocessing'
if resample_only:
	mode = 'Resampling'
table = Table(title='Summary of %s/preprocess.py execution' % app_name,
		box=box.HORIZONTALS,
		show_header=True, header_style='bold magenta')
table.add_column('Mode', justify='center')
table.add_column('# images', justify='center')
table.add_column('Images', justify='center')
table.add_column('Image size', justify='center', no_wrap=True)
table.add_column('Resolution', justify='center')
table.add_row(mode, str(n_imgs),
		str([s1+s2 for s1, s2 in zip(img_fn, img_ext)]),
		str(sz), '%0.4f mm'%img0x.GetSpacing()[0])
console.print(table, justify='center')
console.print('\n')


if resample_only:
	imwrite(img0x, out_path + img_fn[0] + '_x' + img_ext[0])
	rainbow = RainbowHighlighter()
	console.print(rainbow('The first low-res image has been resampled in the high-res lattice'))
	console.print('\n')
	console.print('See it at: [green italic]%s' \
			% out_path + img_fn[0] + '_x' + img_ext[0])
	console.print('\n\n')
	sys.exit(0)

imwrite(img0x, working_path + img_fn[0] + '_x' + img_ext[0])

origin = img0x.GetOrigin()
spacing = img0x.GetSpacing()
direction = img0x.GetDirection()

lr_size = np.zeros([3, n_imgs])
lr_spacing = np.zeros([3, n_imgs])
lr_size[:,0] = np.array(img0.GetSize(), dtype=np.int64)
lr_size[lr_size[:,0]%2!=0,0] -= 1
lr_spacing[:,0] = np.array(img0.GetSpacing())
for ii in track(range(1, n_imgs), '[yellow]Resampling images...'):
	img = imread(working_path + img_fn[ii] + img_ext[ii])
	lr_spacing[:,ii] = np.array(img.GetSpacing())
	lr_size[:,ii] = np.array(img.GetSize())

	lr_size[:,ii] = np.minimum(lr_size[:,ii],
			np.around(spacing / lr_spacing[:,ii] * sz)).astype(np.int64)
	lr_size[lr_size[:,ii]%2!=0,ii] -= 1

	I = resample_img_like(img, img0x)
	imwrite(I, working_path + img_fn[ii] + '_x' + img_ext[ii])

savemat(working_path + 'geo_property.mat', {'sz': sz, 'origin': origin, \
		'spacing': spacing, 'direction': direction})

prefix = [''] + ['reg_'] * (n_imgs - 1)
sufix = ['_x'] * n_imgs
txt = [''.join(s) for s in zip(*[prefix, img_fn, sufix, img_ext])]
with open(working_path + 'data_fn.txt', 'w') as f:
	for ii, fn in enumerate(txt):
		f.write('%s,%s\n' % (fn, 'h_'+img_fn[ii]+'.mat'))

#print('completed step 1')
#print('\t- resample the images')

# step 2: align up all resampled images
for ii in track(range(1, n_imgs), '[magenta]Aligning images...'):
	cmd = 'crlRigidRegistration -t 2 %s%s_x%s %s%s_x%s \
			%sreg_%s_x%s %stfm2_%s.tfm > /dev/null 2>&1' % \
			(working_path, img_fn[0], img_ext[0], \
			working_path, img_fn[ii], img_ext[ii], \
			working_path, img_fn[ii], img_ext[ii], \
			working_path, img_fn[ii])
	os.system(cmd)
#print('completed step 2')
#print('\t- align up all resampled images')

# step 3: create filters for deconvolution
for ii in track(range(0, n_imgs), '[cyan]Creating filters...'):
	fft_win = 1
	max_factor = -np.inf
	for jj in range(0, 3):
		factor = lr_spacing[jj,ii] / spacing[jj]
		if factor > 1:
			# FWHM in the unit of number of pixel and convert it to sigma
			sigma = factor / 2.355
			filter_len = sz[jj] *2
			gw = signal.windows.gaussian(filter_len, std=sigma)
			gw /= np.sum(gw)
			# put it onto 3D space
			shape = np.ones(3, dtype=np.int64)
			shape[jj] = filter_len
			gw = np.reshape(gw, shape)
			gw = np.roll(gw, -filter_len//2, axis=jj)
			# move it to Fourier domain
			GW = np.abs(fftn(gw, [sz[0]*2, sz[1]*2, sz[2]*2]))

			w1_sz = np.array([sz[0]*2, sz[1]*2, sz[2]*2], dtype=np.int64)
			w1_sz[jj] = lr_size[jj,ii]# // 2
			w0_sz = np.array([sz[0]*2, sz[1]*2, sz[2]*2], dtype=np.int64)
			w0_sz[jj] -= lr_size[jj,ii]*2
			w = np.concatenate([np.ones(w1_sz), \
					np.zeros(w0_sz), np.ones(w1_sz)], axis=jj)
			#w = np.abs(fftn(ifftn(w), [sz[0]*2, sz[1]*2, sz[2]*2]))
			#fft_win *= np.transpose(w * GW + 1j * w * GW, axes=[2,1,0])
			if max_factor < factor:
				#fft_win = np.transpose(w * GW + 1j * w * GW, axes=[2,1,0])
				fft_win = np.transpose(GW, axes=[2,1,0])
				max_factor = factor

	savemat(working_path+'h_'+img_fn[ii]+'.mat', {'fft_win': fft_win})

#print('completed step 3')
#print('\t- create filters for deconvolution')

# step 4: volume fusion
z = sitk.GetArrayFromImage(img0x)
L = np.ones_like(z)
for ii in track(range(1, n_imgs), '[medium_purple]Fusing images...'):
	img = imread(working_path + 'reg_' + img_fn[ii] + '_x' + img_ext[ii])
	a = sitk.GetArrayFromImage(img)
	z += a
	L += (a != 0).astype(np.float32)

z[L!=0] = z[L!=0] / L[L!=0]

img_z = np_to_img(z, img0x)
imwrite(img_z, working_path + 'img_mean' + img_ext[0])

#print('completd step 4')
#print('\t- volume fusion')

# rainbow = RainbowHighlighter()
console.print('\n')
console.print('THE PRE-PROCESSING HAS BEEN COMPLETED.')
console.print('\n')
