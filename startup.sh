#!/bin/bash

sudo docker run --rm -v $HOME:/home/app/work/home -p 8888:8888 -p 4005:4005 -it drmeister/cando
