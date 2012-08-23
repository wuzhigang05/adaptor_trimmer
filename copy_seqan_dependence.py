#!/usr/bin/python
import fnmatch
import argparse
import re
import sys
import os
import shutil
import pdb
from tools import line_yielder
parser = argparse.ArgumentParser(
    description= "This program can recursively copy the header files "
    + "that are not in the target directory to the target directory "
    + "from a source directory. This program is written for resolve"
    + "the dependency issue for seqan library used in Adaptor_trimmer"
    + "software. This program has the same function has boost bcp." +
    "Example usage:  python copy_depency.py ~/software/Cplusplus_libs/seqan-trunk/core/include/seqan/ SeqAn1.3/seqan/")

parser.add_argument('source_dir',
    help = "source directory, from which we copy the file to the"
    + "target directory")
parser.add_argument('target_dir', 
    help = "source directory, from which we copy the file to the"
    + "target directory")
parser.add_argument('--pattern', '-p', nargs = '*', default = ["*.h"],
    help = "the file pattern that which we used to glob the files")
args = parser.parse_args()

o = sys.stdout
e = sys.stderr

#counter = 0
def write_list(list, handle = e, sep = "\n"):
  handle.write(sep.join(list) + "\n")

def get_files(dir, pattern_list):
  lst = set()
  for path, dirs, files in os.walk(dir):
    if not path.endswith('/'):
      path += '/'
    files = [path + i for i in files]
    matched_files = [i for i in files for pattern in pattern_list if fnmatch.fnmatch(i, pattern)]
    for i in matched_files:
      lst.add(i)
  return lst

pat_regex = re.compile('#include\s*<\s*(seqan/\S*)\s*>')
def catch_header_from_file(file, c_set):
  for line in line_yielder(file):
    if not line.startswith('#include'):
      continue
    else:
      m = pat_regex.search(line)
      if m:
        c_set.add(m.group(1))
#        if "modifier.h" in line:
#          pdb.set_trace()
def remove_duplicate_parts(dir, file_name):
  paths = dir.rstrip('/').split('/')
  path2 = file_name.split('/')
  for i, j in zip(paths[::-1], path2):
    if i == j:
      paths.remove(i)
    if i != j:
      break
  paths += path2
  str = "/".join(paths) 
#  pdb.set_trace()
  return str 

def check_existence(files):
  existed = []
  not_existed = []
  for i in files:
    if os.path.exists(i):
      existed.append(i)
    else:
      not_existed.append(i)
#  if len(existed) != len(files):
#    e.write("existed: %d not_existed %d\n" % (len(existed), len(not_existed)))
#    e.write("below file was not found:\n")
#    write_list(not_existed)

def get_all_path_valid_headers(source_dir, headers_set):
  s = set()
  for header in headers_set:
    s.add(remove_duplicate_parts(source_dir, header))
  return s


""" This function make sure a directory is created if its not existed 
"""
def check_multiLevel_dir_existense(dir, target_dir):
    paths = os.path.dirname(dir).split('/') # discard the basename
    for i in paths:
      if os.path.exists(os.path.join(target_dir, i)):
        continue
      else:
        os.mkdir(os.path.join(target_dir, i))

def files_not_existed_in_targetdir(file, dir, lst, pat = 'seqan'):
  target_dir = os.path.join(os.getcwd(), dir)
  path_source = file.split('/')
  path_source_sliced = '/'.join(path_source[path_source.index(pat) + 1:])
  valid_path_in_target_file = os.path.join(target_dir, path_source_sliced)
#  pdb.set_trace()
  if not os.path.exists(valid_path_in_target_file):
    lst.add(file)

def copy_file_to_dir(file, dir, pat = 'seqan'):
  target_dir = os.path.join(os.getcwd(), dir)
  path_source = file.split('/')
  path_source_sliced = '/'.join(path_source[path_source.index(pat) + 1:])
  valid_path_in_target_file = os.path.join(target_dir, path_source_sliced)
  not_copied = set()
  if not os.path.exists(valid_path_in_target_file):
    check_multiLevel_dir_existense(path_source_sliced, target_dir)
    try:
      shutil.copy(file, os.path.dirname(valid_path_in_target_file))
      o.write("copied %s to %s\n" % (file, valid_path_in_target_file))
    except:
      e.write("copy file %s failed\n" % (file))
      not_copied.add(file)
      pass
    return not_copied

class Counter:
  i = 0
counter = Counter()

def get_not_existed_files():
  counter.i += 1
  files = get_files(args.target_dir, args.pattern)
  all_headers_set = set()
  for file in files:
    catch_header_from_file(file, all_headers_set)
#  pdb.set_trace()
  all_path_valid_headers = get_all_path_valid_headers(args.source_dir, all_headers_set)
  check_existence(all_path_valid_headers)
  files_not_existed_inTarget = set()
  for file in all_path_valid_headers:
    files_not_existed_in_targetdir(file, args.target_dir, files_not_existed_inTarget)
  e.write("#%d round: below %d files are not in target dir: %s\n" % 
      (counter.i, len(files_not_existed_inTarget), args.target_dir))
  write_list(files_not_existed_inTarget)
#  pdb.set_trace()
  return files_not_existed_inTarget

def process(files_not_existed_inTarget):
  tobedelete = set()
  for file in files_not_existed_inTarget:
    if len(copy_file_to_dir(file, args.target_dir)) == 0:
      tobedelete.add(file)
  for file in tobedelete:
    files_not_existed_inTarget.remove(file)
  new_files_not_existed_inTarget = get_not_existed_files() 
  return new_files_not_existed_inTarget
if __name__ == '__main__':
#  files = get_files(args.target_dir, args.pattern)
#  all_headers_set = set()
#  for file in files:
#    catch_header_from_file(file, all_headers_set)
##  pdb.set_trace()
#  all_path_valid_headers = get_all_path_valid_headers(args.source_dir, all_headers_set)
#  check_existence(all_path_valid_headers)
#  files_not_existed_inTarget = set()
#  for file in all_path_valid_headers:
#    files_not_existed_in_targetdir(file, args.target_dir, files_not_existed_inTarget)
#  e.write("First round: below %d files are not in target dir: %s\n" % 
#      (len(files_not_existed_inTarget), args.target_dir))
#  write_list(files_not_existed_inTarget)
  files_not_existed_inTarget = get_not_existed_files()
  new_files_not_existed_inTarget = process(files_not_existed_inTarget)
  while True:
    if len(new_files_not_existed_inTarget) == 0:
      e.write("All files have been found and copied\n")
      break;
    elif files_not_existed_inTarget == new_files_not_existed_inTarget:
      e.write('*' * 79 + '\n' + "Below file cannot be find in the source dir: %s\n" % (args.source_dir)+
      "Are you sure the source directory is the correct direcotry to use?\n")
      write_list(files_not_existed_inTarget)
      e.write('*' * 79 + '\n')
      break
    

    files_not_existed_inTarget = new_files_not_existed_inTarget
    new_files_not_existed_inTarget = process(new_files_not_existed_inTarget)
