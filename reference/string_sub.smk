"path/{outputname}.ext"

"command {wildcards.outputname} {output}"

expand("path/{filename}.ext", filename=range(10))

expand("{sample,[A-Za-z0-9]+}")
