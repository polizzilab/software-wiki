
# Add the path to our shared installation of Combs so that python knows
# where to import combs2 from

export PYTHONPATH="/n/data1/hms/bcmp/polizzi/lab/programs/design/Combs2/:$PATH"


# We have a shared miniconda executable in our lab shared directory; this block below
# will point your `conda` and your `python` to use it


# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/n/data1/hms/bcmp/polizzi/lab/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/n/data1/hms/bcmp/polizzi/lab/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/n/data1/hms/bcmp/polizzi/lab/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/n/data1/hms/bcmp/polizzi/lab/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
