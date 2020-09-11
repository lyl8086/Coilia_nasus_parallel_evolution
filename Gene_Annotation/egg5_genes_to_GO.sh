#!/bin/bash

bioawk -t '$1!~/^#/ && $7!="" {print $1,$7}' $1
