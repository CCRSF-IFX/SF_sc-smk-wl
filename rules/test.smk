rule touch:
    output: "test.out"
    shell:
        "touch {output}"

rule test:
    default_target: True
    input:
        "test.out"

