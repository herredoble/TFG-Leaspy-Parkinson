import os
import shutil
import distutils.dir_util, distutils.file_util
from pathlib import Path

import nbformat
from nbconvert.exporters import NotebookExporter, export
from nbconvert.preprocessors import ExecutePreprocessor

import jupytext


VERBOSE = True

# +
# Based on your modifications in sources, you can temporary exclude some notebooks to process
# In any case, files at root of sources/ will be copied into challenges/ and solutions/ folders

notebooks_to_process = [
    "sources/TP1_LMM/TP1_LMM.md",
]
#    "sources/TP2_leaspy_beginner/TP2_leaspy_beginner.md",
#    "sources/TP3_advanced_leaspy/TP3_advanced_leaspy.md",

# -

def solution_process(nb, exec_folder):
    """
    Wrapper method that will be called on the solution notebooks, that just runs
    the content of the provided notebook.

    Args:
        nb (nbconvert.NotebookNode): the notebook that will be executed
        exec_folder (Path, or str): working folder where the notebook will be executed.
          This is useful since some notebook involve saving external files.

    Returns:
        A nbconvert.NotebookNode that contains the executed version of the provided nb.
    """
    ep = ExecutePreprocessor(timeout=3600, kernel_name='python3')
    ep.preprocess(nb, {'metadata': {'path': exec_folder}})

    return nb

subversion_config = {
    "solution": {
        "save_folder": lambda x: x.replace("sources", "solutions"),
        "exclude_tags": ["challenge"],
        "extra_post_process": solution_process  # we will execute the solution notebook
    },
    "challenge": {
        "save_folder": lambda x: x.replace("sources", "challenges"),
        "exclude_tags": ["solution"],
    }
}

def process(notebook_path: str):

    # get root folder of notebook and notebook name
    notebook_dirpath = os.path.dirname(notebook_path)
    notebook_name, _ = os.path.splitext(os.path.basename(notebook_path))

    for version, config in subversion_config.items():
        target_folder = config["save_folder"](notebook_dirpath)

        # remove everything that was already present in target folder
        if os.path.isdir(target_folder):
            shutil.rmtree(target_folder)
            os.makedirs(target_folder)

        print(f"Copying {notebook_dirpath}/ into {target_folder}/")
        distutils.dir_util.copy_tree(notebook_dirpath, target_folder, verbose=VERBOSE)

        # delete copied files we do not want to keep (raw source notebooks)
        for p in Path(target_folder).glob(f"{notebook_name}.*"):
            print(f"Deleting {p}")
            p.unlink()

        # delete "hidden" folders (e.g. data/_generation/)
        for p in Path(target_folder).rglob('_*'):
            if p.is_dir():
                print(f"Deleting {p}")
                shutil.rmtree(p)

        # delete potential unwanted folders from release folders
        potential_log_folder = os.path.join(target_folder, "logs")
        potential_output_folder = os.path.join(target_folder, "outputs")
        folders_to_delete = [potential_log_folder, potential_output_folder]
        for folder in folders_to_delete:
            if os.path.exists(folder):
                print(f"Deleting {folder}")
                shutil.rmtree(folder)

        # open original notebook
        original_notebook = jupytext.read(notebook_path)

        new_notebook_path = os.path.join(target_folder, f"{notebook_name}.ipynb")
        print(f"Generating the {version} notebook in {new_notebook_path}")

        exclude_tags_set = set(config["exclude_tags"])
        filtered_notebook = original_notebook.copy()

        # we want to include all cells which are not tagged with an excluded tag
        # so we compute the intersection between each tag and the current set of
        # filtered tags and look whether this intersection is empty to decide
        # note that if `exclude_tags`` contains multiple tags,
        # NONE of them should be present in the cell tags
        filtered_notebook['cells'] = [
            cell for cell in filtered_notebook['cells']
            if exclude_tags_set.intersection(cell['metadata'].get("tags", "default")) == set()
        ]

        # for the solution notebook, we want to execute them, so that we can also
        # save and see the results. This is useful when using Gitlab's notebook
        # viewer, so that we can look at the solutions without having to re-run
        # the notebook
        if config.get("extra_post_process", None) is not None:
            filtered_notebook = config.get("extra_post_process")(filtered_notebook, target_folder)

        exporter = NotebookExporter()
        notebook_string, ressources = export(exporter, nbformat.from_dict(filtered_notebook))

        with open(new_notebook_path, "w+", encoding="utf-8") as notebook_file_handler:
            notebook_file_handler.write(notebook_string)

if __name__ == "__main__":

    # Notebooks to copy & process
    for notebook_path in notebooks_to_process:
        process(notebook_path)

    # Files at root of sources/ are copied into challenges/ and solutions/ folders
    for p in Path('sources').glob('*'):
        if p.is_file() and not p.name.startswith('.'):
            for config in subversion_config.values():
                target_p: str = config["save_folder"](str(p))
                distutils.file_util.copy_file(str(p), target_p, verbose=VERBOSE)




