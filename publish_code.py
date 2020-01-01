import glob
import os

root_dir = os.path.dirname(os.path.realpath(__file__))
all_files = glob.glob(os.path.join(root_dir, "src/*"))

md_head = """
# rosetta-carpet
The [Talbot carpet]() in different languages.
"""

md_template = """
# {LANGUAGE}
```{LANGUAGE}
{CODE}
```

![{LANGUAGE}](carpets/{LANGUAGE}.png)
"""

with open(os.path.join(root_dir, "docs", "index.md"), "w+") as f:
    f.write(md_head)
    for language_file in all_files:
        language_name =  os.path.basename(language_file).rsplit(".")[0]
        print(f"Found language: {language_name}")

        if os.path.isdir(language_file):
            src_files = glob.glob(os.path.join(language_file, "src/*"))
            if len(src_files) == 1:
                language_file = src_files[0]
            else:
                print(f"More than one source file for {language_file}!")

        text = md_template.replace("{LANGUAGE}", language_name)
        with open(language_file, "r") as langf:
            language_code = langf.read()
        text = text.replace("{CODE}", language_code)

        f.write(text)

