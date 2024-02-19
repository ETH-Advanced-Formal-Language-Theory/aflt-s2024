#!/bin/bash

wget https://github.com/Advanced-Formal-Language-Theory/aflt-f2023/releases/download/1.2/pickles.zip
unzip pickles.zip -d autograding_tests
rm pickles.zip

sudo -H pip install -e .
