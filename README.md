# snomed_graph
 
This is a simple Python library for working with [SNOMED CT](https://www.snomed.org).

It has been designed to help newcomers to the terminology to easily load and query it.  The library does not support ECL (Expression Constraint Language), nor should it be considered to be a replacement for an Ontology Server.  Instead, the primary focus of the library is to support analytics and machine learning use cases by making it simple to extract sets of concepts and navigate the relationships between concepts.

The library has three dependencies: 

- `networkx` - any version >= 3.0 should work.
- `pandas`.
- `tqdm`.

Illustrations of the supported functions can be found in the `usage.ipynb` notebook.
