#!/bin/bash
# ===========================================
#  Auto Build & Run Script for TBB + C++17
#  Supports macOS, CentOS7, Ubuntu
# ===========================================

set -e

# --- Detect OS type ---
OS=$(uname -s)
echo "🔍 Detected OS: $OS"

# --- Default compiler ---
CXX=$(which g++)
if [[ -z "$CXX" ]]; then
    echo " g++ compiler not found!"
    echo "➡️  Please install it first:"
    echo "   Ubuntu/Debian: sudo apt install g++"
    echo "   CentOS/RHEL:   sudo yum install gcc-c++"
    exit 1
fi

STD="-std=c++17"
OPT="-O3"

# --- Default TBB paths ---
INCLUDE_PATH=""
LIB_PATH=""

# --- OS-specific configuration ---
if [[ "$OS" == "Darwin" ]]; then
    # macOS (Homebrew)
    if [[ -d "/opt/homebrew/include" ]]; then
        INCLUDE_PATH="-I/opt/homebrew/include"
        LIB_PATH="-L/opt/homebrew/lib"
        echo "Using Homebrew TBB path"
    else
        echo "⚠️  Homebrew TBB not found. Run: brew install tbb"
        exit 1
    fi
else
    # Ubuntu or others
    if [[ -d "/usr/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/include"
        LIB_PATH="-L/usr/lib/x86_64-linux-gnu"
        echo "Using Ubuntu TBB path (/usr/include, /usr/lib/x86_64-linux-gnu)"
    elif [[ -d "/usr/local/include/tbb" ]]; then
        INCLUDE_PATH="-I/usr/local/include"
        LIB_PATH="-L/usr/local/lib"
        echo "Using /usr/local TBB path"
    else
        echo "⚠️  TBB not found. Try: sudo apt install libtbb-dev"
        exit 1
    fi
fi

# --- Source files ---
SRC="Multiple_main_github.cpp functions_github.cpp initial_github.cpp"

# --- Output binary ---
OUT="main"

# --- Build ---
echo "⚙️  Building project..."
$CXX $STD $OPT $INCLUDE_PATH $LIB_PATH $SRC -ltbb -o $OUT

echo "Build completed successfully."

# --- Run program ---
echo "Running program..."
./$OUT
