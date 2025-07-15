#!/bin/bash
cmd="verkko --paths $1 --assembly $2 --hifi "$3" --nano $4 -d $5"
echo $cmd
eval $cmd
