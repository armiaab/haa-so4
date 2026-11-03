"""Editable PhotoCraft shape artwork and matching vector publication exports.

Run with the paper's plotting Python:
  python src/gambar_photocraft.py --import-existing  # first migration only
  python src/gambar_photocraft.py                  # reopen native sources
  python src/gambar_photocraft.py --figure fig4_ansatze_circuits --redraw-circuits

PhotoCraft authors/retains every shape and returns its canonical Bezier paths.
Its PDF exporter is raster-only. Those returned paths are written to SVG and
converted to PDF by VectorCraft, without rasterizing or tracing the artwork.
Matplotlib is used only to obtain Computer Modern / Helvetica glyph outlines.
"""
import argparse
import json
import math
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

from fontTools.misc.transform import Transform
from fontTools.pens.basePen import BasePen
from fontTools.pens.transformPen import TransformPen
from fontTools.svgLib.path import parse_path
from matplotlib import font_manager, rcParams
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path as TextPathCodes
from matplotlib.textpath import TextPath

ROOT = Path(__file__).resolve().parents[1]
FIGURES = ROOT / "figures"
NATIVE = FIGURES / "photocraft"
REVIEW = ROOT / "build" / "photocraft_review"
PHOTOCRAFT = Path(r"C:\Program Files\PhotoCraft\photocraft-cli.exe")
VECTORCRAFT = Path(r"C:\Program Files\VectorCraft\vectorcraft-cli.exe")
PT = 72 / 25.4
SCALE = 4  # PhotoCraft pixels per publication point (288 ppi).
INK, GRAY, PALE, WINE = "#1E293B", "#64748B", "#F1F5F9", "#881337"
SVG = "http://www.w3.org/2000/svg"
XLINK = "http://www.w3.org/1999/xlink"
SIZES = {
    "fig1_vqe_pipeline": (170, 50),
    "fig4_ansatze_circuits": (170, 150),
    "fig5_weyl_chamber": (120, 75),
    "fig6_tangent_directions": (120, 65),
    "fig7_symmetry_manifold": (170, 60),
}


class PhotoCraft:
    def __init__(self):
        self.process = subprocess.Popen(
            [str(PHOTOCRAFT), "serve", "--automation-read-root", str(ROOT),
             "--automation-write-root", str(ROOT)], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, text=True, encoding="utf-8", cwd=ROOT,
        )
        self.serial = 0

    def call(self, method, params=None):
        self.serial += 1
        request = {"id": self.serial, "method": method, "params": params or {}}
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        reply = self.process.stdout.readline()
        if not reply:
            raise RuntimeError("PhotoCraft stopped: " + self.process.stderr.read())
        response = json.loads(reply)
        if not response.get("ok"):
            raise RuntimeError(f"{method}: {response}")
        return response["result"]

    def command(self, command, params=None):
        return self.call("engine.execute", {"command": command, "params": params or {}})

    def close(self):
        self.process.stdin.close()
        self.process.wait(timeout=30)


def knot(point):
    point = list(point)
    return {"anchor": point, "in": point.copy(), "out": point.copy(), "smooth": False}


def preserve_holes(path):
    # PhotoCraft combines subpaths as shape components. Explicit subtraction
    # preserves glyph counters, which SVG represents through contour winding.
    parts=path["subpaths"]
    polygons=[]
    areas=[]
    for part in parts:
        nodes=part["knots"]
        points=[]
        if part["closed"]:
            for a,b in zip(nodes,nodes[1:]+nodes[:1]):
                for step in range(12):
                    t=step/12
                    points.append([(1-t)**3*a["anchor"][k]+3*(1-t)**2*t*a["out"][k]
                                   +3*(1-t)*t*t*b["in"][k]+t**3*b["anchor"][k] for k in (0,1)])
        polygons.append(TextPathCodes(points) if points else None)
        areas.append(sum(a[0]*b[1]-b[0]*a[1] for a,b in zip(points,points[1:]+points[:1]))/2)
    parents=[]
    for i,part in enumerate(parts):
        candidates=[j for j,poly in enumerate(polygons) if poly is not None
                    and abs(areas[j])>abs(areas[i])+1e-6
                    and poly.contains_point(part["knots"][0]["anchor"])] if areas[i] else []
        parents.append(min(candidates,key=lambda j:abs(areas[j])) if candidates else None)
    depth=[]
    for i,part in enumerate(parts):
        ancestor=parents[i]
        winding=0
        level=0
        while ancestor is not None:
            winding+=1 if areas[ancestor]>0 else -1
            level+=1
            ancestor=parents[ancestor]
        sign=1 if areas[i]>0 else -1
        subtract=(level%2==1) if path.get("fillRule")=="evenodd" else (winding!=0 and winding+sign==0)
        part["op"]="subtract" if subtract else "combine"
        depth.append(level)
    path["subpaths"]=[parts[i] for i in sorted(range(len(parts)),key=lambda i:depth[i])]
    return path


class ShapePen(BasePen):
    """Convert SVG/TrueType curves to PhotoCraft's cubic knot representation."""
    def __init__(self):
        super().__init__(None)
        self.subpaths = []
        self.current = None

    def _moveTo(self, p):
        self.current = {"closed": False, "op": "combine", "knots": [knot(p)]}
        self.subpaths.append(self.current)

    def _lineTo(self, p):
        self.current["knots"].append(knot(p))

    def _curveToOne(self, a, b, c):
        self.current["knots"][-1]["out"] = list(a)
        end = knot(c)
        end["in"] = list(b)
        self.current["knots"].append(end)

    def _qCurveToOne(self, control, end):
        start = self.current["knots"][-1]["anchor"]
        a = [start[i] + 2 / 3 * (control[i] - start[i]) for i in range(2)]
        b = [end[i] + 2 / 3 * (control[i] - end[i]) for i in range(2)]
        self._curveToOne(a, b, end)

    def _closePath(self):
        nodes = self.current["knots"]
        if len(nodes) > 1 and math.dist(nodes[0]["anchor"], nodes[-1]["anchor"]) < 1e-6:
            nodes[0]["in"] = nodes.pop()["in"]
        self.current["closed"] = True

    def _endPath(self):
        pass

    def path(self, rule="nonzero"):
        return {"subpaths": self.subpaths, "fillRule": rule, "inverted": False}


def transform_attribute(value):
    result = Transform()
    for name, raw in re.findall(r"(\w+)\s*\(([^)]*)\)", value or ""):
        n = [float(v) for v in re.findall(r"[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?", raw)]
        if name == "translate":
            result = result.translate(n[0], n[1] if len(n) > 1 else 0)
        elif name == "scale":
            result = result.scale(n[0], n[1] if len(n) > 1 else n[0])
        elif name == "matrix":
            result = result.transform(Transform(*n))
        elif name == "rotate":
            cx, cy = n[1:] if len(n) == 3 else (0, 0)
            result = result.translate(cx, cy).rotate(math.radians(n[0])).translate(-cx, -cy)
        else:
            raise ValueError(f"Unsupported SVG transform: {name}")
    return result


def read_svg_shapes(path):
    """Migrate outlined SVG artwork to native shapes, preserving cubic geometry."""
    root = ET.parse(path).getroot()
    ids = {e.attrib["id"]: e for e in root.iter() if "id" in e.attrib}
    shapes = []

    def visit(element, matrix, inherited, label="Artwork"):
        tag = element.tag.rsplit("}", 1)[-1]
        if tag in ("defs", "metadata", "title", "style", "clipPath"):
            return
        style = inherited.copy()
        for key in ("fill", "fill-opacity", "fill-rule", "opacity", "stroke", "stroke-width",
                    "stroke-opacity", "stroke-dasharray", "stroke-linecap", "stroke-linejoin"):
            if key in element.attrib:
                style[key]=element.attrib[key]
        style.update(dict(part.split(":", 1) for part in element.get("style", "").replace(" ", "").split(";") if ":" in part))
        matrix = matrix.transform(transform_attribute(element.get("transform")))
        label = element.get("id", label)
        if tag == "use":
            referenced = ids[element.get(f"{{{XLINK}}}href", element.get("href", ""))[1:]]
            visit(referenced, matrix.translate(float(element.get("x", 0)), float(element.get("y", 0))), style, label)
        elif tag == "path":
            if not element.get("d", "").strip():
                return
            pen = ShapePen()
            parse_path(element.attrib["d"], TransformPen(pen, matrix))
            fill = style.get("fill", "#000000")
            fill = None if fill == "none" else fill
            alpha = float(style.get("fill-opacity", style.get("opacity", 1)))
            if fill and alpha < 1:
                # PhotoCraft accepts normalized RGBA components.
                fill = [int(fill[i:i + 2], 16) / 255 for i in (1, 3, 5)] + [alpha]
            stroke = None
            if style.get("stroke", "none") != "none":
                width = float(style.get("stroke-width", 1)) * SCALE
                stroke = {
                    "width": width, "color": style["stroke"], "align": "center",
                    "cap": style.get("stroke-linecap", "butt"),
                    "join": style.get("stroke-linejoin", "round"),
                    "opacity": 100 * float(style.get("stroke-opacity", style.get("opacity", 1))),
                }
                if style.get("stroke-dasharray", "none") != "none":
                    stroke["dashes"] = [float(v) * SCALE / width for v in re.split(r"[, ]+", style["stroke-dasharray"])]
            item = {"kind": "path", "path": pen.path(style.get("fill-rule", "nonzero")),
                    "fill": fill, "stroke": stroke, "name": label}
            # Consecutive equal styles share one layer (e.g. a line of glyphs).
            if shapes and all(shapes[-1][k] == item[k] for k in ("fill", "stroke")):
                shapes[-1]["path"]["subpaths"].extend(item["path"]["subpaths"])
            else:
                shapes.append(item)
        else:
            for child in element:
                visit(child, matrix, style, label)

    visit(root, Transform().scale(SCALE), {})
    return shapes


class CircuitArtwork:
    def __init__(self):
        self.shapes = []
        self.height = 150
        regular = subprocess.check_output(["kpsewhich", "texgyreheros-regular.otf"], text=True).strip()
        bold = subprocess.check_output(["kpsewhich", "texgyreheros-bold.otf"], text=True).strip()
        font_manager.fontManager.addfont(regular)
        font_manager.fontManager.addfont(bold)
        rcParams["mathtext.fontset"] = "cm"
        self.fonts = {False: FontProperties(fname=regular, math_fontfamily="cm"),
                      True: FontProperties(fname=bold, math_fontfamily="cm")}

    def xy(self, x, y):
        return [x * PT * SCALE, (self.height - y) * PT * SCALE]

    def add(self, kind, name, fill=None, stroke=INK, lw=1.2, dashed=False, **geometry):
        if fill == "white":
            fill = "#FFFFFF"
        edge = None if stroke is None else {"color": stroke, "width": lw * SCALE, "align": "center"}
        if dashed and edge:
            edge["dashes"] = [3 / lw, 2 / lw]
        self.shapes.append({"kind": kind, "name": name, "fill": fill, "stroke": edge, **geometry})

    def line(self, x1, y1, x2, y2, color=INK, lw=1.2, dashed=False):
        # A path stroke, not a filled line shape, retains its exact point width.
        p = {"subpaths": [{"closed": False, "op": "combine",
                            "knots": [knot(self.xy(x1, y1)), knot(self.xy(x2, y2))]}], "fillRule": "nonzero"}
        self.add("path", "Wire", stroke=color, lw=lw, dashed=dashed, path=p)

    def rect(self, x, y, w, h, fill="white", edge=INK, lw=1.2, dashed=False, name="Gate"):
        self.add("rect", name, fill=fill, stroke=edge, lw=lw, dashed=dashed,
                 rect=self.xy(x, y + h) + [w * PT * SCALE, h * PT * SCALE])

    def ellipse(self, x, y, r, fill="white"):
        self.add("ellipse", "CNOT port", fill=fill,
                 rect=self.xy(x - r, y + r) + [2 * r * PT * SCALE] * 2)

    def label(self, x, y, label, size=8, color=INK, bold=False, ha="center"):
        glyphs = TextPath((0, 0), label, size=size, prop=self.fonts[bold])
        box = glyphs.get_extents()
        left = x * PT - (box.x0 + box.width / 2 if ha == "center" else box.x0 if ha == "left" else box.x1)
        baseline = (self.height - y) * PT + (box.y0 + box.height / 2)
        transform = Transform(SCALE, 0, 0, -SCALE, left * SCALE, baseline * SCALE)
        pen = ShapePen()
        target = TransformPen(pen, transform)
        for coordinates, code in glyphs.iter_segments(curves=True):
            pts = list(zip(coordinates[::2], coordinates[1::2]))
            if code == TextPathCodes.MOVETO:
                target.moveTo(pts[0])
            elif code == TextPathCodes.LINETO:
                target.lineTo(pts[0])
            elif code == TextPathCodes.CURVE3:
                target.qCurveTo(*pts)
            elif code == TextPathCodes.CURVE4:
                target.curveTo(*pts)
            elif code == TextPathCodes.CLOSEPOLY:
                target.closePath()
        self.add("path", "Label: " + label, fill=color, stroke=None, path=pen.path())

    def gate(self, x, y, label, w=6, h=5.5, complex_gate=False, edge=INK, size=8):
        self.rect(x-w/2, y-h/2, w, h, PALE if complex_gate else "#FFFFFF",
                  GRAY if complex_gate else edge, .8 if complex_gate else 1.2, complex_gate)
        self.label(x, y, label, size)

    def register(self, left, right, ys, prep=False, trace=False):
        for i, y in enumerate(ys):
            self.line(left, y, right, y)
            self.label(left-2, y, "$a$" if i == 5 else rf"$q_{i}$", ha="right")
        if prep:
            self.gate(left+3, ys[2], "$X$", w=4.8)
        if trace:
            self.gate(right-3.3, ys[-1], r"$\mathrm{Tr}_a$", w=6.6, h=6, size=7.3)

    def repeated(self, left, right, top, bottom):
        self.rect(left, bottom, right-left, top-bottom, None, GRAY, .8, True, "Repeated layer")
        self.label((left+right)/2, top+3.5, r"$\times L$", 9)

    def pair(self, x, top, bottom, star=False):
        self.rect(x-5, bottom-2.5, 10, top-bottom+5, edge=WINE, name="Two-qubit block G")
        if star:
            # Numbered input ports disambiguate nonadjacent wires in a spanning box.
            self.label(x-3.3, top, "0", 5.5, GRAY)
            self.label(x-3.3, bottom, "1", 5.5, GRAY)
        self.label(x+.7 if star else x, (top+bottom)/2, "$G$", 10, WINE)

    def build(self):
        self.rect(0, 0, 170, 150, "#FFFFFF", None, name="Paper")
        for x, y, title in ((2,145,"(a)  HEA"), (88,145,"(b)  HAA"),
                            (2,65,"(c)  Brick block layer"), (88,65,"(d)  Star block layer")):
            self.label(x, y, title, 9, bold=True, ha="left")
        self.label(42,138.5,r"Secondary $R_YR_Z$ variant; primary baseline: $U_3$",7.4,GRAY)
        self.label(128,138.5,"Common structure for CAN and SO(4)",7.4,GRAY)
        ys=[123,116,109,102,95,88]
        self.register(9,81,ys[:5],prep=True)
        self.repeated(17,62,130,90)
        for y in ys[:5]:
            for x in (22,68): self.gate(x,y,"$R_Y$")
            for x in (30,76): self.gate(x,y,"$R_Z$",complex_gate=True)
        for i,x in enumerate((39,45.5,52,58.5)):
            self.line(x,ys[i],x,ys[i+1])
            self.ellipse(x,ys[i],.65,INK)
            self.ellipse(x,ys[i+1],1.65)
            self.line(x-1.65,ys[i+1],x+1.65,ys[i+1])
            self.line(x,ys[i+1]-1.65,x,ys[i+1]+1.65)
        self.label(42,83,r"$X$ on $q_2$ prepares $|\Phi_{\mathrm{HF}}\rangle=|00100\rangle$.",7.5,GRAY)
        self.label(42,77.5,"Final rotations are outside the repeated layer.",7.4,GRAY)

        self.register(97,168,ys,prep=True,trace=True)
        self.repeated(116,151,130,83)
        for y in ys[:5]:
            for x in (109,157): self.gate(x,y,"$V$")
        for y in ys: self.gate(123,y,"$W$")
        self.gate(140,(ys[0]+ys[-1])/2,r"$\mathcal{B}_{\mathcal{P}}$",w=14,h=ys[0]-ys[-1]+5.5,size=13)
        self.label(128,77.5,r"CAN: $V=U_3,\ W=I$; SO(4): $V=I,\ W=R_Y$",7.8)
        self.label(128,71.5,"End rotations have independent parameters.",7.2,GRAY)

        self.label(42,58.5,r"$n_a=0$; two parallel half-layers",8,GRAY)
        self.label(128,58.5,r"$n_a=1$; shared ancilla starts in $|0\rangle$",8,GRAY)
        ys=[51,44,37,30,23,16]
        self.register(9,81,ys[:5])
        self.register(97,168,ys)
        for i,j in ((0,1),(2,3)): self.pair(32,ys[i],ys[j])
        for i,j in ((1,2),(3,4)): self.pair(59,ys[i],ys[j])
        for i,x in enumerate((105,118,131,144,157)):
            self.pair(x,ys[i],ys[-1],star=True)
        self.label(32,16,"Half-layer 1",7.5,GRAY)
        self.label(59,16,"Half-layer 2",7.5,GRAY)
        self.label(128,10.5,"Only ports 0 and 1 participate in each block.",7.3,GRAY)
        self.label(85,4.5,r"$G=\mathrm{CAN}$ or the SO(4) block $B$; panels (c,d) expand $\mathcal{B}_{\mathcal{P}}$.",8)
        return self.shapes


def path_data(path, scale):
    def point(p):
        return " ".join(f"{v / scale:.6f}" for v in p)
    data=[]
    for part in path["subpaths"]:
        nodes=part["knots"]
        if not nodes: continue
        data.append("M "+point(nodes[0]["anchor"]))
        pairs=list(zip(nodes,nodes[1:]))
        if part["closed"]: pairs.append((nodes[-1],nodes[0]))
        for a,b in pairs:
            if a["out"] == a["anchor"] and b["in"] == b["anchor"]:
                data.append("L "+point(b["anchor"]))
            else:
                data.append("C "+" ".join(point(p) for p in (a["out"],b["in"],b["anchor"])))
        if part["closed"]: data.append("Z")
    return " ".join(data)


def fill_attributes(fill, prefix):
    if fill is None: return {prefix:"none"}
    if isinstance(fill,str): return {prefix:fill}
    if isinstance(fill,list):
        return {prefix:"#"+"".join(f"{round(v*255):02x}" for v in fill[:3]),
                prefix+"-opacity":str(fill[3] if len(fill)>3 else 1)}
    raise ValueError(f"Unexpected PhotoCraft paint: {fill}")


def export_shapes(shapes,name,size):
    ET.register_namespace("",SVG)
    root=ET.Element(f"{{{SVG}}}svg",{
        "width":f"{size[0]}mm","height":f"{size[1]}mm",
        "viewBox":f"0 0 {size[0]*PT:.8f} {size[1]*PT:.8f}","version":"1.1",
    })
    ET.SubElement(root,f"{{{SVG}}}title").text=name
    ET.SubElement(root,f"{{{SVG}}}desc").text="PhotoCraft shape layers; canonical Bezier paths exported without rasterization."
    for i,shape in enumerate(shapes):
        attrs={"id":f"shape-{i}","d":path_data(shape["path"],SCALE),
               "fill-rule":shape["path"]["fillRule"],**fill_attributes(shape["fill"],"fill")}
        stroke=shape.get("stroke")
        if stroke:
            attrs.update(fill_attributes(stroke["fill"],"stroke"))
            attrs.update({"stroke-width":str(stroke["width"]/SCALE),
                          "stroke-linecap":stroke["cap"],"stroke-linejoin":stroke["join"],
                          "stroke-opacity":str(stroke["opacity"]/100)})
            if stroke.get("dashes"):
                attrs["stroke-dasharray"]=" ".join(str(v*stroke["width"]/SCALE) for v in stroke["dashes"])
        p=ET.SubElement(root,f"{{{SVG}}}path",attrs)
        ET.SubElement(p,f"{{{SVG}}}title").text=shape["name"]
    svg=FIGURES/(name+".svg")
    ET.ElementTree(root).write(svg,encoding="utf-8",xml_declaration=True)
    result=subprocess.run([str(VECTORCRAFT),"convert",str(svg),str(FIGURES/(name+".pdf"))],
                          check=True,capture_output=True,text=True,encoding="utf-8")
    report=json.loads(result.stdout)
    (REVIEW/(name+"_export.json")).write_text(json.dumps(report,indent=2),encoding="utf-8")
    if report.get("warnings") or report.get("importWarnings"):
        raise RuntimeError(f"Vector export warnings for {name}: {report}")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--import-existing",action="store_true")
    parser.add_argument("--svg-directory",type=Path,default=FIGURES)
    parser.add_argument("--redraw-circuits",action="store_true",help="Rebuild Figure 4 from its gate definitions")
    parser.add_argument("--redesign",action="store_true",help="Rebuild Figures 1, 5, 6, and 7 from the concept designs")
    parser.add_argument("--figure",choices=SIZES)
    args=parser.parse_args()
    NATIVE.mkdir(exist_ok=True)
    REVIEW.mkdir(parents=True,exist_ok=True)
    for name,size in SIZES.items():
        if args.figure and name!=args.figure: continue
        app=PhotoCraft()
        try:
            native=NATIVE/(name+".pcraft")
            if args.redesign:
                from desain_teori import BUILDERS
            redesign=args.redesign and name in BUILDERS
            redraw=redesign or name=="fig4_ansatze_circuits" and (args.redraw_circuits or not native.exists())
            if redraw or args.import_existing:
                shapes=BUILDERS[name]() if redesign else CircuitArtwork().build() if name=="fig4_ansatze_circuits" else read_svg_shapes(args.svg_directory/(name+".svg"))
                app.call("doc.new",{"width":round(size[0]*PT*SCALE),"height":round(size[1]*PT*SCALE),
                                    "resolution":72*SCALE,"background":"transparent","name":name})
                canonical=[]
                for i,shape in enumerate(shapes):
                    if shape["kind"]=="path" and shape.get("fill") is not None:
                        shape["path"]=preserve_holes(shape["path"])
                    canonical.append(app.command("shape.create",shape))
                    if (i+1)%50==0: print(f"{name}: {i+1}/{len(shapes)} PhotoCraft layers",flush=True)
                app.call("doc.save",{"path":native.relative_to(ROOT).as_posix()})
            else:
                app.call("doc.open",{"path":native.relative_to(ROOT).as_posix()})
                document=app.call("doc.inspect")
                canonical=[app.command("shape.info",{"layer":layer["id"]})
                           for layer in reversed(document["layers"]) if layer["kind"]=="Shape"]
            app.call("doc.save",{"path":(REVIEW/(name+".png")).relative_to(ROOT).as_posix()})
            export_shapes(canonical,name,size)
            print(f"{name}: {len(canonical)} editable PhotoCraft shapes; PDF/SVG exported",flush=True)
        finally:
            app.close()


if __name__=="__main__":
    main()
