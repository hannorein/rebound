import rebound
import inspect
import docstring_to_markdown
def convert_code_blocks(doc):
    new_doc = ""
    lines = doc.split("\n")
    first = True
    for line in lines:
        if first:
            if line[:3]==">>>":
                first = False
                new_doc += "```python\n"
                new_doc += line[3:]+"\n"
            else:
                new_doc += line+"\n"
        else:
            if line[:3]==">>>":
                new_doc += line[3:]+"\n"
            else:
                new_doc += "```\n"
                new_doc += line+"\n"
                first = True
    if first==False:
        new_doc += "```\n"

    return new_doc

def render_class(cls, functions=None):
    d = "## Class `"+cls+"`\n"
    d += convert_code_blocks(inspect.cleandoc(eval(cls).__doc__))
    for function in functions:
        f = getattr(eval(cls),function)
        d += "## Function `"+cls+"."+function+"`\n"
        d += convert_code_blocks(inspect.cleandoc(f.__doc__))

    return d

print(render_class("rebound.Simulation",["copy"]))

