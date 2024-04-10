import sys

def main(arg1, arg2):
    editHtmlFile(arg1, arg2)

def editHtmlFile(filename, outputFile):
    f = open(filename)
    lines = f.readlines()
    f.close()
    f = open(outputFile, 'w')
    for line in lines:
        if 'data:application/pdf' in line:
            line = line.replace('<img src=', '<embed src=')
            line = line.replace('/>\n', ' width="600" height="625"/>\n')
        f.write(line)
    f.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
