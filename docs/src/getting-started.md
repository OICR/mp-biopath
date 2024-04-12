Getting Started
This chapter covers everything you need to know to start using MP-BioPath to analyze your files. This includes installing MP-BioPath and running a few example commands.

Installation
To simplify the installation process, a Docker container has been generated. All you need to do is pull the latest version of the container.

If you're not familiar with Docker, it provides an environment with all the dependencies pre-installed to run the tool.

Docker
Docker is a program that runs on your computer. You'll need to install it on your OS X, Windows, or Linux system.

Instructions for installing Docker can be found here: Docker Installation Guide

Docker Image
To obtain the latest Docker container, use the following command:

bash
Copy code:
docker pull oicr/mp-biopath:1.0.x

After running this command, check that you have the container by typing:

bash
Copy code:
docker images

You should see the image "oicr/mp-biopath" in your list of images.

Running the Container
To run a command line from within the container, use:

bash
Copy code:
make run-bash

This command places you inside the container. Note that the filesystem inside the container is separate from the host filesystem. The "-v" flag maps a directory outside the container to a folder inside the container. In this example, it maps the current directory to a "data" directory inside the container. Changes to files made within this directory will not be lost when exiting the container.

To run commands inside the container from the host, use:

bash
Copy code:
docker run mpbiopath julia test/run.jl

This command runs the tests for the tool. For running custom commands, refer to the "Running MP-BioPath" section.

Custom Setup
To run the tool on your own data, you need to create evidence files in the specified format. You also need to obtain the Reactome pathways in pairwise interaction format from Reactome Pathway Analysis.

After obtaining the required files, create a custom configuration file. MP-BioPath will use this configuration file to run and produce results in the specified folder.

Refer to specific sections of this document for details on developing files in the correct format.
