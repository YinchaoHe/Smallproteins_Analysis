from typing import List
from pafpy import PafFile

def main():
    cazyfamilies = []
    path = "output/minimap2.paf"
    with open(path) as fileobj:
        paf = PafFile(fileobj)
        for record in paf:
            cazyfamily = record.qname.split('/')[0]
            if cazyfamily in cazyfamilies:
                continue
            else:
                cazyfamilies.append(cazyfamily)
    print(len(cazyfamilies))


if __name__ == '__main__':
    main()
