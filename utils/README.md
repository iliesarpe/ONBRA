## How to remap a dataset

- `python remap.py <filename.txt> <SAVEIDS>`
- `<filename.txt>` is the path to the file in format `src dst tim` to be remapped to conform to the input format
- `<SAVEIDS>` is either `0` then the IDS of the original nodes will not be stored or `1` that stores the map of the original node IDs

The script takes an arbitrary temporal network and remaps it to be used for input to the `onbra` executable. As an example if you execute the script with `python remap.py sample.txt 1` (you find `sample.txt` in the current folder) you will find in this same folder the file `sample-remapped.txt` which will correspond to the remapped dataset and `sample-IDMAP.txt` which will be a map between the original IDs and the novel IDs assigned by the script to the nodes
