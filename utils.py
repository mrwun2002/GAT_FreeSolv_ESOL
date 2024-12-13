# quick fix
# from https://github.com/rdkit/rdkit/blob/master/rdkit/Chem/Draw/__init__.py

import math
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw.canvasbase import CanvasBase
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem.Draw.MolDrawing import Font


class MPLCanvas(CanvasBase):
    def __init__(self, fig, axes, font):
        self.size = (200, 200)
        self._dpi = 200
        self._figure = fig
        self._axes = axes
        self.points = []
        self.font = Font(**font)

    def rescalePt(self, p1):
        return [float(p1[0]) / self._dpi, float(self.size[1] - p1[1]) / self._dpi]

    def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
        canvas = self._axes
        p1 = self.rescalePt(p1)
        p2 = self.rescalePt(p2)
        if color2 and color2 != color:
            mp = (p1[0] + p2[0]) / 2.0, (p1[1] + p2[1]) / 2.0
            canvas.add_line(
                Line2D((p1[0], mp[0]), (p1[1], mp[1]), color=color, **kwargs)
            )
            canvas.add_line(
                Line2D((mp[0], p2[0]), (mp[1], p2[1]), color=color2, **kwargs)
            )
        else:
            canvas.add_line(
                Line2D((p1[0], p2[0]), (p1[1], p2[1]), color=color, **kwargs)
            )
        self.points.append(p1)
        self.points.append(p2)

    def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
        import re

        pos = self.rescalePt(pos)

        ax = self._axes
        fig = self._figure

        text = re.sub(r"<sub>(.*?)</sub>", r"$_{\1}$", text)
        text = re.sub(r"<sup>(.*?)</sup>", r"$^{\1}$", text)
        text = re.sub(r"<.*?>", "", text)

        orientation = kwargs.get("orientation", "E")
        halign = "center"
        valign = "center"
        x, y = pos

        if orientation != "C":
            c = ax.text(
                x,
                y,
                text[0],
                verticalalignment="baseline",
                horizontalalignment=halign,
                weight=font.weight,
                size=self.font.size,
                family=self.font.face,
            )
            bb = c.get_window_extent(fig.canvas.get_renderer())
            (x0, y0), (x1, y1) = ax.transData.inverted().transform(bb)
            w = x1 - x0
            h = y1 - pos[1]
            c.remove()

            alpha = 0.0

            if orientation == "E":
                halign = "left"
                valign = "center"
                x = pos[0] - w / 2.0 + alpha * w
                y = pos[1] #- (pos[1] - y0) / 2.0
            elif orientation == "W":
                halign = "right"
                valign = "center"
                x = pos[0] + w / 2.0 - alpha * w
                y = pos[1] #- (pos[1] - y0) / 2.0
            elif orientation == "N":
                valign = "bottom"
                x = pos[0]
                y = y0 + alpha * h
            elif orientation == "S":
                valign = "top"
                x = pos[0]
                y = pos[1] + (pos[1] - y0) / 2.0 - alpha * h

        annot = ax.text(
            x,
            y,
            text,
            color=color,
            verticalalignment=valign,
            horizontalalignment=halign,
            weight=font.weight,
            size=self.font.size,
            family=self.font.face,
            bbox=dict(boxstyle="square,pad=0", fc="w", lw=0),
        )

        fig.canvas.draw()
        bb = annot.get_window_extent()
        bb2 = bb.transformed(self._axes.transData.inverted())
        self.points.append([bb2.x0, bb2.y0])
        self.points.append([bb2.x1, bb2.y1])

        try:
            bb = annot.get_window_extent()
            w, h = bb.width, bb.height
            tw, th = ax.transData.inverted().transform((w, h))
        except AttributeError:
            tw, th = 0.1, 0.1  # <- kludge
        return (tw, th, 0)

    def addCanvasPolygon(self, ps, color=(0, 0, 0), **kwargs):
        canvas = self._axes
        ps = [self.rescalePt(x) for x in ps]
        for x in ps:
            self.points.append(x)
        canvas.add_patch(Polygon(ps, linewidth=0, facecolor=color))

    def addCanvasDashedWedge(
        self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs
    ):
        canvas = self._axes
        dash = (3, 3)
        pts1 = self._getLinePoints(p1, p2, dash)
        pts2 = self._getLinePoints(p1, p3, dash)
        pts1 = [self.rescalePt(p) for p in pts1]
        pts2 = [self.rescalePt(p) for p in pts2]

        if len(pts2) < len(pts1):
            pts2, pts1 = pts1, pts2
        for i in range(len(pts1)):
            if color2 and color2 != color:
                mp = (pts1[i][0] + pts2[i][0]) / 2.0, (pts1[i][1] + pts2[i][1]) / 2.0
                canvas.add_line(
                    Line2D(
                        (pts1[i][0], mp[0]), (pts1[i][1], mp[1]), color=color, **kwargs
                    )
                )
                canvas.add_line(
                    Line2D(
                        (mp[0], pts2[i][0]), (mp[1], pts2[i][1]), color=color2, **kwargs
                    )
                )
            else:
                canvas.add_line(
                    Line2D(
                        (pts1[i][0], pts2[i][0]),
                        (pts1[i][1], pts2[i][1]),
                        color=color,
                        **kwargs
                    )
                )
            self.points.append(pts1[i])
            self.points.append(pts2[i])


def DrawMolToMPL(
    mol,
    fig,
    ax,
    font={"face": "sans", "size": 12},
    style="golden",
    kekulize=True,
    wedgeBonds=True,
    imageType=None,
    fitImage=False,
    options=None,
    **kwargs
):
    if not mol:
        raise ValueError("Null molecule provided")

    size = fig.get_size_inches() * fig.dpi
    canvas = MPLCanvas(fig, ax, font)
    if options is None:
        options = DrawingOptions()
        options.bgColor = None
    if fitImage:
        options.dotsPerAngstrom = int(min(size) / 10)
    options.wedgeDashedBonds = wedgeBonds
    drawer = MolDrawing(canvas=canvas, drawingOptions=options)
    omol = mol
    if kekulize:
        mol = Chem.Mol(mol.ToBinary())
        Chem.Kekulize(mol)

    if not mol.GetNumConformers():
        AllChem.Compute2DCoords(mol)

    ax.plot([0, 1], [0, 1], linestyle="None")
    drawer.AddMol(mol, **kwargs)
    omol._atomPs = drawer.atomPs[mol]
    for k, v in omol._atomPs.items():
        omol._atomPs[k] = canvas.rescalePt(v)

    x_max = max([x for x, y in canvas.points])
    x_min = min([x for x, y in canvas.points])
    y_max = max([y for x, y in canvas.points])
    y_min = min([y for x, y in canvas.points])
    # [x_min, x_min, x_max, x_max, x_min], [y_min, y_max, y_max, y_min, y_min]

    if style == "tight":
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
    elif style == "golden":
        gr = (1 + math.sqrt(5)) / 2
        w = x_max - x_min
        h = y_max - y_min
        if w > gr * h:
            ycent = y_min + h / 2.0
            yh = w / gr
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(ycent - yh / 2, ycent + yh / 2)
        else:
            xcent = x_min + w / 2.0
            xh = h * gr
            ax.set_xlim(xcent - xh / 2, xcent + xh / 2)
            ax.set_ylim(y_min, y_max)
    ax.set_aspect("equal")
    ax.set_axis_off()
    return (x_min, x_max), (y_min, y_max)


def draw_smiles(smiles_string, fig = None, ax = None):
    molecule = Chem.MolFromSmiles(smiles_string)

    if fig == None:
        fig = plt.gcf()
    if ax == None:
        ax = fig.gca()
    
    DrawMolToMPL(molecule, fig, ax)



def my_draw_networkx_edge_labels(
        #https://stackoverflow.com/questions/22785849/drawing-multiple-edges-between-two-nodes-with-networkx
    G,
    pos,
    edge_labels=None,
    label_pos=0.5,
    font_size=10,
    font_color="k",
    font_family="sans-serif",
    font_weight="normal",
    alpha=None,
    bbox=None,
    horizontalalignment="center",
    verticalalignment="center",
    ax=None,
    rotate=True,
    clip_on=True,
    rad=0
):
    """Draw edge labels.

    Parameters
    ----------
    G : graph
        A networkx graph

    pos : dictionary
        A dictionary with nodes as keys and positions as values.
        Positions should be sequences of length 2.

    edge_labels : dictionary (default={})
        Edge labels in a dictionary of labels keyed by edge two-tuple.
        Only labels for the keys in the dictionary are drawn.

    label_pos : float (default=0.5)
        Position of edge label along edge (0=head, 0.5=center, 1=tail)

    font_size : int (default=10)
        Font size for text labels

    font_color : string (default='k' black)
        Font color string

    font_weight : string (default='normal')
        Font weight

    font_family : string (default='sans-serif')
        Font family

    alpha : float or None (default=None)
        The text transparency

    bbox : Matplotlib bbox, optional
        Specify text box properties (e.g. shape, color etc.) for edge labels.
        Default is {boxstyle='round', ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0)}.

    horizontalalignment : string (default='center')
        Horizontal alignment {'center', 'right', 'left'}

    verticalalignment : string (default='center')
        Vertical alignment {'center', 'top', 'bottom', 'baseline', 'center_baseline'}

    ax : Matplotlib Axes object, optional
        Draw the graph in the specified Matplotlib axes.

    rotate : bool (deafult=True)
        Rotate edge labels to lie parallel to edges

    clip_on : bool (default=True)
        Turn on clipping of edge labels at axis boundaries

    Returns
    -------
    dict
        `dict` of labels keyed by edge

    Examples
    --------
    >>> G = nx.dodecahedral_graph()
    >>> edge_labels = nx.draw_networkx_edge_labels(G, pos=nx.spring_layout(G))

    Also see the NetworkX drawing examples at
    https://networkx.org/documentation/latest/auto_examples/index.html

    See Also
    --------
    draw
    draw_networkx
    draw_networkx_nodes
    draw_networkx_edges
    draw_networkx_labels
    """
    import matplotlib.pyplot as plt
    import numpy as np

    if ax is None:
        ax = plt.gca()
    if edge_labels is None:
        labels = {(u, v): d for u, v, d in G.edges(data=True)}
    else:
        labels = edge_labels
    text_items = {}
    for (n1, n2), label in labels.items():
        (x1, y1) = pos[n1]
        (x2, y2) = pos[n2]
        (x, y) = (
            x1 * label_pos + x2 * (1.0 - label_pos),
            y1 * label_pos + y2 * (1.0 - label_pos),
        )
        pos_1 = ax.transData.transform(np.array(pos[n1]))
        pos_2 = ax.transData.transform(np.array(pos[n2]))
        linear_mid = 0.5*pos_1 + 0.5*pos_2
        d_pos = pos_2 - pos_1
        rotation_matrix = np.array([(0,1), (-1,0)])
        ctrl_1 = linear_mid + rad*rotation_matrix@d_pos
        ctrl_mid_1 = 0.5*pos_1 + 0.5*ctrl_1
        ctrl_mid_2 = 0.5*pos_2 + 0.5*ctrl_1
        bezier_mid = 0.5*ctrl_mid_1 + 0.5*ctrl_mid_2
        (x, y) = ax.transData.inverted().transform(bezier_mid)

        if rotate:
            # in degrees
            angle = np.arctan2(y2 - y1, x2 - x1) / (2.0 * np.pi) * 360
            # make label orientation "right-side-up"
            if angle > 90:
                angle -= 180
            if angle < -90:
                angle += 180
            # transform data coordinate angle to screen coordinate angle
            xy = np.array((x, y))
            trans_angle = ax.transData.transform_angles(
                np.array((angle,)), xy.reshape((1, 2))
            )[0]
        else:
            trans_angle = 0.0
        # use default box of white with white border
        if bbox is None:
            bbox = dict(boxstyle="round", ec=(1.0, 1.0, 1.0), fc=(1.0, 1.0, 1.0))
        if not isinstance(label, str):
            label = str(label)  # this makes "1" and 1 labeled the same

        t = ax.text(
            x,
            y,
            label,
            size=font_size,
            color=font_color,
            family=font_family,
            weight=font_weight,
            alpha=alpha,
            horizontalalignment=horizontalalignment,
            verticalalignment=verticalalignment,
            rotation=trans_angle,
            transform=ax.transData,
            bbox=bbox,
            zorder=1,
            clip_on=clip_on,
        )
        text_items[(n1, n2)] = t

    ax.tick_params(
        axis="both",
        which="both",
        bottom=False,
        left=False,
        labelbottom=False,
        labelleft=False,
    )

    return text_items