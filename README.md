<img src="https://github.com/Katzlab/PhyloToL-6/blob/main/Other/Katzlab.png">

> :warning: ** This is currently being dockerised! ** Use the container at your own peril!

# **EukPhylo version 1.0** 
EukPhylo is an updated version of the PhyloToL pipeline from the [Katz Lab](https://www.science.smith.edu/katz-lab/) at Smith College. EukPhylo is a phylogenomic toolkit for processing transcriptomic and genomic data that includes novel phylogeny-informed contamination removal techniques. See our [Wiki](https://github.com/Katzlab/EukPhylo/wiki) for more information on installation and usage!


## Dockerfile

The [docker file](https://github.com/Katzlab/EukPhylo/blob/Docker/PTL1/Dockerfile.txt) for part 1 can be executed with:

```bash
# Build the container
docker build -f Dockerfile.txt . --tag eukphylo


# Current command is:
docker run -it \
    --mount type=bind,src=$(pwd)/databases,dst=/Databases \
    --mount type=bind,src=$(pwd)/input_data,dst=/Input_data \
    --mount type=bind,src=$(pwd)/output_data,dst=/Output_data \
    eukphylo
```

An example for running the dockerfile that takes in an OGlist, taxonlist, and R2Gs as input. It also requires an Output folder.
> :warning: Do not change the "dst=/$(path)", only change "src=$(pwd)"
```
docker run -it \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/OG_list.txt,dst=/EukPhylo/PTL2listofOGs.txt \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/taxon_list.txt,dst=/EukPhylo/PTL2taxon_list.txt \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/R2G,dst=/Input_data \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/Output_data,dst=/Output_data \
eukphylo
```



The [docker file](https://github.com/Katzlab/EukPhylo/blob/Docker/PTL2/Dockerfile.txt) for part 2 can be executed with:

```bash
# Build the container
docker build -f Dockerfile.txt . --tag eukphylo


# Current command is:
docker run -it \
    --mount type=bind,src=$(pwd)/databases,dst=/Databases \
    --mount type=bind,src=$(pwd)/input_data,dst=/Input_data \
    --mount type=bind,src=$(pwd)/output_data,dst=/Output_data \
    eukphylo
```

An example for running the dockerfile that takes in an OGlist, taxonlist, and R2Gs as input. It also requires an Output folder.
> :warning: Do not change the "dst=/$(path)", only change "src=$(pwd)"
```
docker run -it \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/OG_list.txt,dst=/EukPhylo/PTL2listofOGs.txt \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/taxon_list.txt,dst=/EukPhylo/PTL2taxon_list.txt \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/R2G,dst=/Input_data \
--mount type=bind,src=/Users/gani/phylotol_ms/Docker/PT2/Output_data,dst=/Output_data \
eukphylo
```


After development, GitHub CICD workflows can be added to automatically build and release the dockerfile for the end user.
