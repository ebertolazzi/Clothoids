#!/bin/bash

astyle --style=1tbs --indent=spaces=2 --attach-namespaces --attach-classes \
       --indent-namespaces --indent-preproc-block --indent-preproc-define \
       --indent-col1-comments --pad-oper --pad-comma --pad-paren \
       --align-pointer=middle --align-reference=middle \
       --keep-one-line-blocks --keep-one-line-statements \
       --convert-tabs $1
