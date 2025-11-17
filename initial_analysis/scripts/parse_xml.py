import xml.etree.ElementTree as ET

def parse_xml(filename):
    # Load and parse the XML file
    tree = ET.parse(filename)
    root = tree.getroot()
    
    # Function to recursively print the XML structure
    def print_elements(element, indent=0):
        print(' ' * indent + f'<{element.tag}>: {element.text.strip()}')
        for child in element:
            print_elements(child, indent + 2)
    
    # Print the root and its child elements
    print_elements(root)

if __name__ == "__main__":
    filename = "hmdb_metabolites.xml"  # replace with your XML file path
    parse_xml(filename)
