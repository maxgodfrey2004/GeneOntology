.PHONY: clean install

install:
	@echo "> Downloading go-basic.obo ..."
	curl -o data/go-basic.obo http://current.geneontology.org/ontology/go-basic.obo
	@echo "> Downloading goa_human.gaf.gz ..."
	curl -o data/goa_human.gaf.gz http://current.geneontology.org/annotations/goa_human.gaf.gz
	@echo "> Unzipping goa_human.gaf ..."
	gzip -d data/goa_human.gaf.gz

clean:
	rm data/*
