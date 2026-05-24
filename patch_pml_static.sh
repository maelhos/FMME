#!/bin/bash

find . -type f \( -name "Makefile" -o -name "makefile" \) -exec sed -i 's/-lntl/\/usr\/local\/lib\/libntl.a/g' {} +