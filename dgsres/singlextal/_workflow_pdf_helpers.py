import pylatex

def initReportDoc(fname):
    geometry_options = {"right": "2cm", "left": "2cm"}
    doc = pylatex.Document(fname, geometry_options=geometry_options)
    return doc

def verbatim(s):
    import pylatex
    return pylatex.utils.NoEscape( r'\begin{verbatim}' + s + '\end{verbatim}' )

def slice_info_section(sl, doc):
    import pylatex
    with doc.create(pylatex.Section('Slice info')):
        lines = [
            'Name: %s' % sl.name,
            'hkl0: %s' % (sl.hkl0,),
            'projection: %s' % (sl.hkl_projection,),
            ]
        s = '\n'.join(lines)
        doc.append(verbatim(s))
    return
