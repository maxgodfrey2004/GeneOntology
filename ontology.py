from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from collections import namedtuple

GOEntry = namedtuple('GOEntry', 'goid name def_ is_a')
HumanAnnotation = namedtuple('HumanAnnotation', 'symbol relation goid')

# Stores a basic version of the gene ontology. Keys are GOids.
GO_DATA = {}

# Contains Human Genetic Annotations
HUMAN_DATA = []

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
        in_term = True
        goid = None
        name = None
        def_ = None
        is_a = []
      elif line == '':
        entry = GOEntry(goid, name, def_, is_a)
        in_term = False
        GO_DATA[entry.goid] = entry
      elif in_term:
        if line.startswith('is_a'):
          is_a.append(line.split()[1])
        elif line.startswith('def'):
          def_ = line.split('"')[1]
        elif line.startswith('id'):
          goid = line.split()[1]
        elif line.startswith('name'):
          name = ' '.join(line.split()[1:])

def load_goa_human_gaf():
  '''Loads the information in data/goa_human.gaf
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

print('Loading go-basic.obo')
load_go_obo()
print('Loading goa_human.gaf')
load_goa_human_gaf()

def get_cellular_component(symbol):
  '''Returns the GOids of the cellular components a gene pertains to.
  '''
  ids = set()
  for annotation in HUMAN_DATA:
    if annotation.symbol == symbol and annotation.relation == 'located_in':
      ids.add(annotation.goid)
  return [(goid, GO_DATA[goid].def_) for goid in ids]

def get_biological_process(symbol):
  '''Returns the GOids of the biological processes a gene pertains to.
  '''
  ids = set()
  for annotation in HUMAN_DATA:
    if annotation.symbol == symbol and annotation.relation == 'involved_in':
      ids.add(annotation.goid)
  return [(goid, GO_DATA[goid].def_) for goid in ids]

def get_molecular_function(symbol):
  '''Returns the GOids of the molecular functions a gene pertains to.
  '''
  ids = set()
  for annotation in HUMAN_DATA:
    if annotation.symbol == symbol and annotation.relation == 'enables':
      ids.add(annotation.goid)
  return [(goid, GO_DATA[goid].def_) for goid in ids]
  

if __name__ == '__main__':
  print(get_biological_process('PCNA'))
  print(get_cellular_component('PCNA'))
  print(get_molecular_function('PCNA'))