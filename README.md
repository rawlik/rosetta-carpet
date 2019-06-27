# rosetta-carpet

The Talbot carpet in many languages!

The source code and resulting carpets are hosted on the webpage [https://rawlik.github.io/rosetta-carpet/](https://rawlik.github.io/rosetta-carpet/).

The script `publish_code.py` automatically generates the website based on the source code in the `src` directory. It assumes that the code is a single file named after the programming language: `LANGUAGE_NAME.extension`. It also expects to find a `docs/carpets/LANGUAGE_NAME.png` file. The code is printed on the website using the template:
````markdown
```LANGUAGE_NAME
content of the src/LANGUAGE_NAME.extension file
```
````