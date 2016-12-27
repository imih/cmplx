#!/bin/bash
#
# Downloads and installs the latest nightly snapshot of igraph into
# the home directory of the user

###########################################################################

set -e

###########################################################################

readlink_that_works()
{
    DIR="$1"
    READLINK=`which readlink`
    if [ x$READLINE != x ]; then
        readlink -f "$DIR"
    else
        # Mac OS X does not have readlink -f. Meh.
        perl -e 'use Cwd "abs_path";print abs_path(shift)' "$DIR"
    fi
}

###########################################################################

identify_distro()
{
    echo -n "identifying operating system... "
    OS=`uname -s`
    echo "$OS"
    if [ "x$OS" = "xLinux" ]; then
        echo -n "identifying distribution... "
        if [ -f /etc/debian_version ]; then
            echo "Debian or Ubuntu"
            DISTRO="Debian"
        else
            echo "unknown"
            DISTRO="generic"
        fi

        INSTALL_DEPS_FUNC="${OS}_${DISTRO}"
    else
        INSTALL_DEPS_FUNC="generic"
    fi

    INSTALL_DEPS_FUNC="install_dependencies_${INSTALL_DEPS_FUNC}"
}

###########################################################################

install_dependencies()
{
    echo "1. Installing dependencies"
    echo "=========================="
    echo ""

    identify_distro
    echo ""

    ${INSTALL_DEPS_FUNC}
    echo ""
}

install_dependencies_Linux_Debian()
{
    DEBIAN_DEPS="gcc bzr flex bison automake autoconf libxml2-dev \
                 libglpk-dev libgmp3-dev python-dev"

    REQUIRED=`dpkg -l ${DEBIAN_DEPS} 2>/dev/null | grep -v "^.i" | sed '1,5d; s/^....\([^ ]*\).*/\1/g' | tr '\n' ' '`
    REQUIRED2=`dpkg -l ${DEBIAN_DEPS} 2>&1 >/dev/null | grep "No packages found" | sed -e 's/No packages found matching \([^ ]*\)\..*/\1/g' | tr '\n' ' '`
    REQUIRED="${REQUIRED} ${REQUIRED2}"

    if [ "x${REQUIRED}" = x -o "x${REQUIRED}" = "x " ]; then
        echo "All required packages are installed."
        return
    fi

    if [ ${IS_ROOT} = 1 ]; then
        echo "Installing required packages..."
        apt-get install ${REQUIRED}
    else
        echo "You are not running this script with root privileges and the"
        echo "following required packages are not installed:"
        for REQ in ${REQUIRED}; do
            echo "  - ${REQ}"
        done

        while [ x$RESPONSE != xy ]; do
            echo ""
            echo "Do you want me to try and install them? You will need root"
            echo "privileges to perform this task and you will be asked for"
            echo "the root password. [Y/n]"
            read RESPONSE
            if [ x$RESPONSE = xy -o x$RESPONSE = xY ]; then
                RESPONSE=y
            elif [ x$RESPONSE = xn -o x$RESPONSE = xN ]; then
                RESPONSE=n
                echo "Continuing without installing missing dependencies."
                return
            else
                echo "Please respond with either y or n."
            fi
        done

        sudo apt-get install ${REQUIRED}
    fi
}

install_dependencies_Linux_generic()
{
    echo "/!\\ Could not identify exact Linux distribution."
    echo "    Let's hope that all the dependencies are installed."
}

install_dependencies_generic()
{
    echo "/!\\ Could not identify the operating system."
    echo "    Let's hope that all the dependencies are installed."
}

###########################################################################

check_dependencies()
{
    echo "2. Checking dependencies"
    echo "========================"
    echo ""

    MISSING=
    for DEPS in ${DEPS_REQUIRED_TOOLS}; do
        OK=0
        for DEP in `echo ${DEPS}|tr '|' ' '`; do
            echo -n "checking ${DEP}... "
            if which $DEP; then
                OK=1
                break
            else
                echo "not found!"
            fi
        done
        if [ $OK = 0 ]; then
            MISSING="${MISSING} ${DEPS}"
        fi
    done
    if [ "x$MISSING" != x ]; then
        echo ""
        echo "/!\\ The following required dependencies are missing:"
        for DEPS in ${MISSING}; do
            DEPS=`echo ${DEPS}|sed -e 's/|/ or /g'`
            echo "    - ${DEPS}"
        done
        exit 1
    fi

    # Compare the prefixes of autoconf, automake and libtool
    AC_PREFIX=`dirname \`which autoconf\``
    AM_PREFIX=`dirname \`which automake\``
    LT_PREFIX=`dirname \`which libtool\``
    if [ "x${AC_PREFIX}" != "x${AM_PREFIX}" ]; then
        echo ""
        echo "/!\\ WARNING: the prefixes of autoconf and automake are"
        echo "    different. You may run into problems during installation."
        echo "    Proceed with caution."
    fi
    if [ "x${AC_PREFIX}" != "x${LT_PREFIX}" ]; then
        echo ""
        echo "/!\\ WARNING: the prefixes of autoconf and libtool are"
        echo "    different. You may run into problems during installation."
        echo "    Proceed with caution."
    fi

    # Optional dependencies
    echo -n "checking Python... "
    PYTHON=`which python`
    if [ "x${PYTHON}" = x ]; then
        echo "not found, skipping Python interface"
        HAVE_PYTHON=0
    else
        echo "${PYTHON}"
        HAVE_PYTHON=1
    fi

    # Determining Python version
    if [ ${HAVE_PYTHON} = 1 ]; then
        echo -n "checking Python version... "
        PYVERSION=`${PYTHON} -c "import platform; print('.'.join(platform.python_version_tuple()[:2]))"`
        echo "${PYVERSION}"

        echo -n "checking Python development headers... "
        PYHEADER="`${PYTHON} -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_inc())"`/Python.h"
        if [ -f "${PYHEADER}" ]; then
            echo "${PYHEADER}"
        else
            echo "not found, skipping Python interface"
            HAVE_PYTHON=0
        fi
    fi

    echo ""
}

###########################################################################

checkout_source()
{
    echo "3. Checking out igraph from Launchpad"
    echo "====================================="
    echo ""
    echo "Repository URL: ${REPO}"
    echo "Please wait while the source code is being checked out..."
    mkdir -p "${SRCDIR}/igraph"
    cd "${SRCDIR}/igraph"
    rm -rf ${BRANCH}
    bzr co ${REPO} ${BRANCH}
    echo ""
}

###########################################################################

compile_core()
{
    echo "4. Compiling C core of igraph"
    echo "============================="
    echo ""
    cd "${SRCDIR}/igraph/${BRANCH}"
    ./bootstrap.sh
    mkdir build
    cd build
    ../configure --prefix="${PREFIX}" $@
    make
    echo ""
}

###########################################################################

install_core()
{
    echo "5. Installing C core of igraph"
    echo "=============================="
    echo ""
    cd "${SRCDIR}/igraph/${BRANCH}/build"
    make install
    echo ""
}

###########################################################################

compile_python_interface()
{
    echo "6. Compiling Python interface"
    echo "============================="
    echo ""
    cd "${SRCDIR}/igraph/${BRANCH}/interfaces/python"

    # Patching setup.py to disable pkg-config and fall back to setup.pcfg
    sed -i -e 's/pkg-config igraph/pkg-config noigraphplease/g' setup.py

    # Build the extension
    ${PYTHON} setup.py build
    echo ""
}

###########################################################################

install_python_interface()
{
    echo "7. Installing Python interface"
    echo "=============================="
    echo ""
    cd "${SRCDIR}/igraph/${BRANCH}/interfaces/python"

    # Create the site-packages dir in the prefix if it does not exist
    SITE_PKG_DIR="${PREFIX}/lib/python${PYVERSION}/site-packages"
    if [ ! -d "${SITE_PKG_DIR}" ]; then
        mkdir -p ${SITE_PKG_DIR}
    fi

    # Ensure that SITE_PKG_DIR is in PYTHONPATH and install the extension
    PYTHONPATH="${SITE_PKG_DIR}:${PYTHONPATH}" \
        ${PYTHON} setup.py install --prefix=${PREFIX}
    echo ""

    # Overwrite the launch script
    if [ $OS = Darwin ]; then
        LD_LIBRARY_PATH_VARNAME="DYLD_LIBRARY_PATH"
    else
        LD_LIBRARY_PATH_VARNAME="LD_LIBRARY_PATH"
    fi
    cat >${PREFIX}/bin/igraph <<EOF
#!/bin/bash
if [ "x\${PYTHONPATH}" = x ]; then
	export PYTHONPATH=${SITE_PKG_DIR}
else
	export PYTHONPATH=${SITE_PKG_DIR}:\${PYTHONPATH}
fi

if [ "x\${${LD_LIBRARY_PATH_VARNAME}}" = x ]; then
	export ${LD_LIBRARY_PATH_VARNAME}=${PREFIX}/lib
else
	export ${LD_LIBRARY_PATH_VARNAME}=${PREFIX}/lib:\${${LD_LIBRARY_PATH_VARNAME}}
fi

python -m igraph.app.shell
EOF
}

###########################################################################

usage()
{
    echo "Usage: $0 [options]"
    echo ""
    echo "Checks out the latest snapshot of igraph from Launchpad, compiles"
    echo "it and tries to install it using a given prefix."
    echo ""
    echo "The default invocation without any arguments will check out the"
    echo "development branch and compile and install it into the home"
    echo "directory of the user."
    echo ""
    echo "Options:"
    echo ""
    echo "    -b BRANCH   specifies the branch to be checked out. Default: trunk"
    echo "    -d DIR      compile igraph in DIR/igraph/BRANCH. Default: ~/src"
    echo "    -h          shows this help message"
    echo "    -p PREFIX   install igraph into PREFIX. Default: ~"
    echo "    -v VERSION  compile and install the given VERSION of igraph."
}

###########################################################################

# Set up defaults
BRANCH="trunk"
SRCDIR=""
PREFIX=""

###########################################################################

# Parse command line options

while getopts ":b:d:hp:v:" opt; do
    case $opt in
    b)
        BRANCH="$OPTARG"
        ;;
    d)
        SRCDIR=$(readlink_that_works "$OPTARG")
        ;;
    h)
        usage
        exit 0
        ;;
    p)
        PREFIX=$(readlink_that_works "$OPTARG")
        ;;
    v)
        BRANCH="$OPTARG-main"
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        exit 1
        ;;
    :)
        echo "Option -$OPTARG requires an argument." >&2
        exit 1
        ;;
    esac
done

shift $((OPTIND-1))

###########################################################################

# Set up dependent variables

DEPS_REQUIRED_TOOLS="bzr gcc g++ make automake autoconf libtoolize|glibtoolize flex bison"
REPO="https+urllib://code.launchpad.net/igraph/${BRANCH}/"

if [ "x${SRCDIR}" = x ]; then
    SRCDIR="${HOME}/src"
fi
if [ "x${PREFIX}" = x ]; then
    PREFIX="${HOME}"
fi

if [ `id -u` = 0 ]; then
    IS_ROOT=1
else
    IS_ROOT=0
fi

###########################################################################

echo "Checking out branch: ${BRANCH}"
echo "Compiling igraph in: ${SRCDIR}/igraph/${BRANCH}"
echo "Installing igraph into: ${PREFIX}"
echo ""

install_dependencies
#check_dependencies

checkout_source
compile_core
install_core

if [ ${HAVE_PYTHON} = 1 ]; then
    compile_python_interface
    install_python_interface
fi
