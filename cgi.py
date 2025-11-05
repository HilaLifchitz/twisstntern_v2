from html import escape as html_escape

def escape(s, quote=True):
    return html_escape(s, quote=quote)
