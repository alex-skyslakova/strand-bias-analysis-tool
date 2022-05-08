INSTALL_JELLYFISH=false

ifeq (, $(shell which jellyfish)) # check jellyfish is not installed
 	INSTALL_JELLYFISH=true
endif

jellyfish-linux:
	@if [ $(INSTALL_JELLYFISH) = "true" ]; then\
		echo "jellyfish not installed. installing...";\
        sudo apt-get update;\
        sudo apt-get install jellyfish;\
    fi

jellyfish-mac:
	@if [ $(INSTALL_JELLYFISH) = "true" ]; then\
		echo "jellyfish not installed. installing...";\
        brew install jellyfish;\
    fi

.PHONY: jellyfish