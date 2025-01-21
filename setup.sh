#!/bin/bash

# 获取当前脚本所在目录
PARP_DIR="$(cd "$(dirname "$0")" && pwd)"

# 打印当前脚本所在目录
echo "PWD: $PARP_DIR"
export PARP_DIR
export PYTHONPATH="$PYTHONPATH:$PARP_DIR"
