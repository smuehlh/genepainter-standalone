#!/usr/bin/env python
# encoding: utf-8
"""
Tool to render a newick phb tree file into svg.

Licence: MPI
Author: Marcel Hellkamp <marc@gsites.de>
"""

import os
import sys
import re

def main(argv):
    from optparse import OptionParser

    parser = OptionParser()
    parser.add_option("--fontsize", default=12.0, type="float", metavar="PX",
                      help="Size of label font")
    parser.add_option("--spacing", default=2.0, type="float", metavar="PX",
                      help="Space between labels")
    parser.add_option("--textwidth", default=100, type="int", metavar="PX",
                      help="Horizontal space to reserve for labels")
    parser.add_option("--maxwidth", default=1024, type="int", metavar="PX",
                      help="Maximum width of entire image")
    parser.add_option("--minstep", default=1, type="float", metavar="PX",
                      help="Minimum width of one branch.")

    (opt, args) = parser.parse_args()

    tree = NewickTree()

    with open(args[0], 'r') as fp:
        tree.parse(fp)

    with open(args[1], 'w') as fp:
        svg = render_svg(tree, fontsize=opt.fontsize, spacing=opt.spacing, textwidth=opt.textwidth, maxwidth=opt.maxwidth, minstep=opt.minstep)
        fp.write(svg)

    return 0


def next(s, search):
    ''' Returns the index of the first occurance of multible search terms in
        a string or -1 on error (no match). Use a list or a multiline string
        to search for multible words. Use a string to search for multible
        chars.
    '''
    return min([x for x in [s.find(x) for x in search] if x >= 0] or [-1])


class Tree(object):
    """ Represents a Tree, a subtree or a leaf. """

    def __init__(self, parent):
        '''Implements a tree'''
        self.label = ''
        self.parent = parent
        self.childs = []
        self._allchilds = []
        self._leafs = []
        self._path = []
        if self.parent:
            self.parent.childs.append(self)

    def left(self):
        """ First (left) child """
        return self.childs[0]

    def right(self):
        ''' Last (right) child '''
        return self.childs[-1]

    def isroot(self):
        return (self.parent == None)

    def isleaf(self):
        return (not self.childs)

    def isempty(self):
        ''' Empty nodes can have a parent, but no childs or label '''
        return not (self.childs or self.label)

    def path(self):
        if not self._path:
            if self.isroot():
                self._path = [self]
            else:
                self._path = self.parent.path() + [self]
        return self._path

    def allchilds(self):
        if not self._allchilds:
            for c in self.childs:
                self._allchilds += c.allchilds() + [c]
        return self._allchilds

    def leafs(self):
        '''Leafs, left first'''
        if not self._leafs:
            for c in self.childs:
                if c.isleaf():
                    self._leafs += [c]
                else:
                    self._leafs += c.leafs()
        return self._leafs

    def add(self, child):
        for c in child.allchilds():
            c._path = []
        for p in self.path():
            self._allchilds = []
            self._leafs = []
        self.childs.append(child)
        child.parent = self

    def remove(self, child):
        for c in child.allchilds():
            c._path = []
        for p in self.path():
            self._allchilds = []
            self._leafs = []
        self.childs.remove(child)
        child.parent = None

    def swap(self, p):
        """Swap this node with a given node"""
        my_parent = self.parent
        his_parent = p.parent
        if my_parent:
            my_parent.remove(self)
        if his_parent:
            his_parent.remove(p)
        if my_parent:
            my_parent.add(p)
        if his_parent:
            his_parent.add(self)

    def insert(self, p):
        """ Insert this node into another tree, replacing the given node and
            adding the removed node as child. """
        if self.parent:
            self.parent.remove(self)

        if p.parent:
            p.parent.add(self)
            p.parent.remove(p)
        self.add(p)

    def reset(self):
        """ Resets metadata up and down the tree. Only touches meta
            information that is affected by the current node """
        for p in self.path():
            self._leafs = []
            self._allchilds = []
        self._path = []
        for c in self.allchilds():
            c._path = []

    def render(self, order=0):
        print ' ' * order, self.label
        for c in self.childs:
            c.render(order + 1)

    def nodes(self):
        yield self
        for n in self.allchilds():
            yield n


class WeightedTree(Tree):
    def __init__(self, parent=None):
        Tree.__init__(self, parent)
        self.weight = 0.0

    def maxdepth(self):
        return max([l.position()[0] for l in self.tree.leafs()])


class NewickTree(WeightedTree):
    def parse(self, io, s=''):
        '''Parses a tree from a newick file and returns the rest '''
        while True:
            while io and len(s) < 128:
                n = io.read(128)
                if n:
                    s += ''.join(n.split())
                else:
                    break
            if s[0] == "(":  # Start of new subtree
                s = s[1:]
                sub = self.__class__(self)  #  Create new subtree
                s = sub.parse(io, s)  #  continue with new subtree
            elif s[0] == ',':  # A comma seperates childs.
                s = s[1:]
                if self.isempty():  #  last node was empty -> ',,' or '(,'
                    return s  #   Close it and continue with parent
                sub = self.__class__(self)  # Start a new subtree
                s = sub.parse(io, s)
            else:  # expect a label label or a :dist
                if s[0] == ')':  # '()x:y' is the same as 'x:y'
                    s = s[1:]
                    if self.isempty():
                        self.parent.childs.remove(self)
                        return s
                x = next(s, ',)')  # Labels end with , (we are a left child) or ')' (we are a last child)
                if x < 0:
                    label = s
                    s = ''
                else:
                    label = s[:x]
                    s = s[x:]
                if label:  # Label found
                    x = next(label, ':')
                    if x == 0:
                        self.weight = abs(float(label[1:]))
                    elif x > 0:
                        self.label = label[:x]
                        self.weight = abs(float(label[x + 1:]))
                    else:
                        self.label = label
                return s
        self.reset()
        return s

    def export(self):
        label = self.label
        if not label:
            label = ''
        if self.weight:
            label += ':%f' % self.weight
        if self.isleaf():
            return label
        else:
            return '(' + ','.join([c.export() for c in self.childs]) + ')' + label




def tree_layout(tree):
    """ Returns three lists of nodes and lines prepared for a renderer:
    node  = (x, y, width, text, isleaf) ordered by y --> A-------B Text
    line  = Vertical line used to connect childs to parent node (x, y, height) ordered by x
    """
    nodes = []
    lines = []
    lc = [0.5]
    
    minweight = min(n.weight for n in tree.allchilds() if n.weight)
    for n in tree.allchilds():
        if not n.weight:
            n.weight = minweight
    
    # print tree.export()

    def go_down(node, xpos):
        # print xpos
        if node.isleaf():
            ''' a---b text'''
            nodes.append((xpos, lc[0], node.weight, node.label, True))
            lc[0] += 1
            return lc[0] - 1
        else:
            ''' a
                |
            c--(d)
                |
                b
            '''
            cpos = [go_down(c, xpos + node.weight) for c in node.childs]
            ax = xpos + node.weight
            ay = min(cpos)
            bx = ax
            by = max(cpos)
            lines.append((ax, ay, by - ay))
            if not node.isroot():
                cx = xpos
                cy = (ay + by) / 2
                dx = ax
                dy = cy
                nodes.append((cx, cy, node.weight, node.label, False))
                return cy

    go_down(tree, 0.0)
    nodes.sort(lambda a, b: cmp(a[1], b[1]))  # Order by y inplace
    lines.sort()  # Order by x inplace
    return nodes, lines


def render_svg(tree, fontsize=12, spacing=2, textwidth=100, 
               maxwidth=sys.maxint, minstep=1):
    nodes, lines = tree_layout(tree)
    minvalue = min(n[2] for n in nodes if n[2] > 0)
    sx, sy = (float(minstep) / minvalue, fontsize + spacing)
    image_width = max(x + w for x, y, w, text, isleaf in nodes) * sx + textwidth + 1 # +1: avoid text cutting of at edge
    image_height = sum(1 for node in nodes if node[4]) * sy + 1 # +1: avoid text cutting of at edge
    
    if image_width > maxwidth:
        sx = sx * (maxwidth - textwidth) / image_width

    svg = ['<svg xmlns="http://www.w3.org/2000/svg"'
           ' xmlns:xlink="http://www.w3.org/1999/xlink"'
           ' width="%d" height="%d">'
           '<g  transform="translate(.5,.5)">' % (image_width, image_height)]

    def node(name, *a, **kwa):
        attr = ' '.join('%s="%s"' % (k.strip('_').replace('_','-'),int(v) if isinstance(v, (int, float)) else v) for (k,v) in kwa.iteritems())
        if a:
            svg.append('<%s %s>' % (name, attr))
            svg.extend(a)
            svg.append('</%s>' % name)
        else:
            svg.append('<%s %s/>' % (name, attr))
    
    smallfont = int(fontsize* 0.85)
    
    linestyle = "stroke:rgb(0,0,0);stroke-width:1"
    leaflabel = "font-size: %dpx"%fontsize
    branchlabel = "font-size: %dpx"%smallfont
    
    for x, y, h in lines:
        node('line', x1=x*sx, y1=y*sy, x2=x*sx, y2=(y+h)*sy, style=linestyle)
    for x, y, w, text, isleaf in nodes:
        node('line', x1=x*sx, y1=y*sy, x2=(x+w)*sx, y2=y*sy, style=linestyle)
        node('circle', cx=(x+w)*sx, cy=y*sy, r=2, fill="black")
        if not text: continue
        m = re.match('(.*)_([0-9\-]+)([^_]+)?_([0-9\-]+)([^_]+)?$', text)
        cl = 'taxon'
        if m:
            text, v1, c1, v2, c2 = m.groups()
            c1 = c1 or 'black'
            c2 = c2 or 'black'
            node('text', v2, x=(x+w)*sx-3, y=y*sy+smallfont-1, text_anchor="end", style=branchlabel, fill=c2)
            node('text', v1, x=(x+w)*sx-3-len(v2)*smallfont, y=y*sy+smallfont-1, text_anchor="end", style=branchlabel, fill=c1)
            if c1 == 'green' and v1 != '0':
                cl = cl+' taxon-intron-gain'
            else:
                cl = 'taxon'   
        if isleaf:
            node('text', text.replace('.', ' '), x=(x+w)*sx + 4, y=y*sy+fontsize/3, style=leaflabel, _class=cl)
        elif text:
            node('text', text, x=(x+w)*sx-3, y=y*sy-3,
                     text_anchor="end", style=branchlabel, _class=cl)

    svg.append('</g></svg>')
    return '\n'.join(svg)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
