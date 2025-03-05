<img src="https://github.com/Katzlab/PhyloToL-6/blob/main/Other/Katzlab.png">

> :warning: ** This is currently being dockerised! ** Use the container at your own peril!

# **EukPhylo version 1.0** 
EukPhylo is an updated version of the PhyloToL pipeline from the [Katz Lab](https://www.science.smith.edu/katz-lab/) at Smith College. EukPhylo is a phylogenomic toolkit for processing transcriptomic and genomic data that includes novel phylogeny-informed contamination removal techniques. See our [Wiki](https://github.com/Katzlab/EukPhylo/wiki) for more information on installation and usage!


## Dockerfile

The docker file can be executed with:

```bash
# Build the container
docker build -f Dockerfile.txt . --tag eukphylo:1

# Get the container IMAGE_ID
docker image list

# Current command is:
docker run -it \
    --mount type=bind,src=$(pwd)/databases,dst=/Databases \
    --mount type=bind,src=$(pwd)/input_data,dst=/Input_data \
    --mount type=bind,src=$(pwd)/output_data,dst=/Output_data \
    eukphylo
```

After development, GitHub CICD workflows can be added to automatically build and release the dockerfile for the end user.
