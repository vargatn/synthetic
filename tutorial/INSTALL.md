# Automatic install

Use the supplied [install script](../install.sh)

Be sure to have anaconda setup on your machine

    conda 23.1.0

run the install script in **interactive mode** this is important for the environment management

    bash -i install.sh

The script aims to set up a new conda environment for you named `synth`,
with all dependencies and the `synthetic` package installed

you can activate it as

    conda activate synth

and deactivate it as

    conda deactivate

If something goes wrong, you can delete it as

    conda env remove -n synth