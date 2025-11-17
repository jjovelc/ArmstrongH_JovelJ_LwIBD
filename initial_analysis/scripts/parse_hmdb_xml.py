import xml.etree.ElementTree as ET

def extract_info(element):
    """
    Extract required fields from an XML element.
    If any field is missing, it returns 'Unknown' for that field.
    """
    fields = [
        "{http://www.hmdb.ca}name",
        "{http://www.hmdb.ca}accession",
        ".//{http://www.hmdb.ca}kingdom",
        ".//{http://www.hmdb.ca}super_class",
        ".//{http://www.hmdb.ca}class",
        ".//{http://www.hmdb.ca}sub_class",
        ".//{http://www.hmdb.ca}molecular_framework"
    ]
    return [element.find(field).text if element.find(field) is not None else 'Unknown' for field in fields]

def main():
    tree = ET.parse("hmdb_metabolites.xml")
    root = tree.getroot()

    with open("results.tsv", "w") as file:
        file.write(f'Metabolite\tAccession\tKingdom\tSuper Class\tClass\tSub Class\tMolecular Framework\n')

        metabolite_elements = list(root.findall(".//{http://www.hmdb.ca}metabolite"))

        # Limit to first 100 entries
        for i, element in enumerate(metabolite_elements):
            if i >= 100:
                break
            
            info = extract_info(element)
            file.write("\t".join(info) + "\n")

if __name__ == "__main__":
    main()
