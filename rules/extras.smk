#ruleorder: gz_compress > gz_extract > faToTwoBit > twoBitToFa
ruleorder: gz_extract > twoBitToFa > faToTwoBit

rule twoBitToFa:
    input: "data/2bit/{genome}.2bit"
    output: "data/genome/{genome}.fa"
    conda: "../envs/conda_twoBitToFa.yaml"
    shell:
        "twoBitToFa {input} {output}"

rule faToTwoBit:
    input: "data/genome/{genome}.fa"
    output: protected("data/2bit/{genome}.2bit")
    conda: "../envs/conda_faToTwoBit.yaml"
    shell:
        "faToTwoBit {input} {output}"

rule gz_extract:
    input: "{file}.fa.gz"
    output: "{file}.fa"
    shell: "gzip -d {input}"

#rule gz_compress:
#    input: "{file}.fa"
#    output: "{file}.fa.gz"
#    wildcard_constraints:
#        extension="[A-Za-z0-9]+(?!gz)"
#    shell: "gzip -k {input}"
