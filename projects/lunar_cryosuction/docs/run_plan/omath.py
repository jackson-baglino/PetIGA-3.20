"""LaTeX -> MathML (latex2mathml) -> OMML (Word's own mathml2omml.xsl)."""
from latex2mathml.converter import convert
from lxml import etree

XSL = "/Applications/Microsoft Word.app/Contents/Resources/mathml2omml.xsl"
_X = etree.XSLT(etree.parse(XSL))


def omml(latex):
    return _X(etree.fromstring(convert(latex))).getroot()
