from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from collections import namedtuple

import sys

GOEntry = namedtuple('GOEntry', 'goid name def_ is_a')
HumanAnnotation = namedtuple('HumanAnnotation', 'symbol relation goid')

# Stores a basic version of the gene ontology. Keys are GOids.
GO_DATA = {}

# Contains Human Genetic Annotations.
HUMAN_DATA = []

# Contains PCNA interactors.
INTERACTORS = []

def load_go_obo():
  '''Loads the information in data/go-basic.obo.
  '''
  global GO_DATA
  in_term = False
  goid = None
  name = None
  def_ = None
  is_a = None
  with open('data/go-basic.obo', 'r') as infile:
    for line in infile:
      line = line.strip()
      if line == '[Term]':
        # A new entity is being defined.
        in_term = True
        goid = None
        name = None
        def_ = None
        is_a = []
      elif line == '':
        # The entity has finished being defined: add it to the GO_DATA map.
        entry = GOEntry(goid, name, def_, is_a)
        in_term = False
        GO_DATA[entry.goid] = entry
      elif in_term:
        # There is more information to parse about the current entity.
        if line.startswith('is_a'):
          is_a.append(line.split()[1])
        elif line.startswith('def'):
          def_ = '"' + line.split('"')[1] + '"'
        elif line.startswith('id'):
          goid = line.split()[1]
        elif line.startswith('name'):
          name = ' '.join(line.split()[1:])

def load_goa_human_gaf():
  '''Loads the information in data/goa_human.gaf.
  '''
  global HUMAN_DATA
  with open('data/goa_human.gaf', 'r') as infile:
    for line in infile:
      line = line.strip()
      if line[0] == '!':
        continue  # Part of the heading
      line = line.split()
      assert line[0] == 'UniProtKB'
      # We are interested in the following indices of `line`:
      #   [2]: Gene symbol
      #   [3]: Relational qualifier
      #   [4]: GOid
      annotation = HumanAnnotation(line[2], line[3], line[4])
      HUMAN_DATA.append(annotation)

def load_interactors():
  '''Loads the interactor information in the static folder.
  '''
  global INTERACTORS
  unique_interactors = set()
  with open('static/reactome-results.txt', 'r') as infile:
    for line in infile:
      line = line.strip()
      unique_interactors.add(line.split()[1])
  with open('static/string-results.txt', 'r') as infile:
    for line in infile:
      line = line.strip()
      unique_interactors.add(line)
  INTERACTORS = list(unique_interactors)

def get_cellular_component(symbol):
  '''Returns the GOids of the cellular components a gene pertains to.
  '''
  ids = set()
  for annotation in HUMAN_DATA:
    if annotation.symbol == symbol and annotation.relation == 'located_in':
      ids.add(annotation.goid)
  return list(ids)

def get_biological_process(symbol):
  '''Returns the GOids of the biological processes a gene pertains to.
  '''
  ids = set()
  for annotation in HUMAN_DATA:
    if annotation.symbol == symbol and annotation.relation == 'involved_in':
      ids.add(annotation.goid)
  return list(ids)

def get_molecular_function(symbol):
  '''Returns the GOids of the molecular functions a gene pertains to.
  '''
  ids = set()
  for annotation in HUMAN_DATA:
    if annotation.symbol == symbol and annotation.relation == 'enables':
      ids.add(annotation.goid)
  return list(ids)

try:
  print('Loading go-basic.obo')
  load_go_obo()
  print('Loading goa_human.gaf')
  load_goa_human_gaf()
except FileNotFoundError:
  print('Error: GO files not found. Run `make install`.', file=sys.stderr)
  sys.exit(1)
try:
  print('Loading interactors')
  load_interactors()
except FileNotFoundError:
  print('Error: Static Interactor files not found. Pull them from github.',
        file=sys.stderr)

if __name__ == '__main__':
  # Sanity checks for annotation information retrieval
  # print(get_biological_process('PCNA'))
  # print(get_cellular_component('PCNA'))
  # print(get_molecular_function('PCNA'))

  # Retrieve all molecular functions in the system
  mf_id_counts = {goid:1 for goid in get_molecular_function('PCNA')}
  for interactor in INTERACTORS:
    for goid in get_molecular_function(interactor):
      mf_id_counts[goid] = mf_id_counts.get(goid, 0) + 1
  entries = [f'{goid},{GO_DATA[goid].def_},{mf_id_counts[goid]}' for goid in mf_id_counts]
  entries.sort(key=lambda x: int(x.split(',')[-1]), reverse=True)
  entries = ['GOid,Molecular Function Annotation Description,Frequency'] + entries
  print(*entries, sep='\n', file=open('results/molecular_function.csv', 'w'))

  # Retrieve all cellular components in the system
  cc_id_counts = {goid:1 for goid in set(get_cellular_component('PCNA'))}
  for interactor in INTERACTORS:
    for goid in get_cellular_component(interactor):
      cc_id_counts[goid] = cc_id_counts.get(goid, 0) + 1
  entries = [f'{goid},{GO_DATA[goid].def_},{cc_id_counts[goid]}' for goid in cc_id_counts]
  entries.sort(key=lambda x: int(x.split(',')[-1]), reverse=True)
  entries = ['GOid,Cellular Component Annotation Description,Frequency'] + entries
  print(*entries, sep='\n', file=open('results/cellular_component.csv', 'w'))

  # Retrieve all biological processes in the system
  bp_id_counts = {goid:1 for goid in get_biological_process('PCNA')}
  for interactor in INTERACTORS:
    for goid in get_biological_process(interactor):
      bp_id_counts[goid] = bp_id_counts.get(goid, 0) + 1
  entries = [f'{goid},{GO_DATA[goid].def_},{bp_id_counts[goid]}' for goid in bp_id_counts]
  entries.sort(key=lambda x: int(x.split(',')[-1]), reverse=True)
  entries = ['GOid,Biological Process Annotation Description,Frequency'] + entries
  print(*entries, sep='\n', file=open('results/biological_process.csv', 'w'))

