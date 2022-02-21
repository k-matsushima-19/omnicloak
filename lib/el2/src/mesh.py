#!/usr/bin/env python3
# coding:utf-8
import math
from math import sin, cos, sqrt
import csv
import numpy as np
import random

# ul = None # 周期長
# x_p0 = None # 周期境界の左端の座標

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __add__(self, other):
        return Point(self.x+other.x, self.y+other.y)

    def __sub__(self, other):
        return Point(self.x-other.x, self.y-other.y)

# Nodeは座標と自身のuniqueな番号を持つ
class Node(Point):
    def __init__(self, x, y):
        Point.__init__(self, x, y)
        self.n = None

# Linkは両端の2つのNodeと自身のuniqueな番号を持つ
class Link:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2
        self.n = None

# Shapeは分割数とNode数，Link数を持つ
class Shape:
    def __init__(self, nn, nl, domain):
        self.nn = nn
        self.nl = nl
        self.nodes = []
        self.links = []
        self.domain = domain

    # LinkにNodeを割り当てる
    def gen_link(self, closed=True):
        # 2次元の場合は順番につなげればいい
        self.links = [Link(self.nodes[i],self.nodes[i+1]) for i in range(self.nn-1)]
        if closed: self.links.append(Link(self.nodes[self.nn-1],self.nodes[0]))

# CanvasはShapeと現在の累計Node数，Link数を持つ
class Canvas:
    def __init__(self):
        self.shapes = []
        self.nn = 0 # Nodeの総数
        self.nl = 0 # Linkの総数
        self.inners = []
        self.nip = 0 # 内点の総数

    def add(self, shape):
        # shapeが持つすべてのNodeとLinkにuniqueな番号を振る
        tmp = 0
        for node in shape.nodes:
            node.n = self.nn + tmp
            tmp = tmp + 1
        self.nn = self.nn + shape.nn
        tmp = 0
        for link in shape.links:
            link.n = self.nl + tmp
            tmp = tmp + 1
        self.nl = self.nl + shape.nl
        self.shapes.append(shape)

    def add_inner(self, inner):
        self.inners.append(inner)
        self.nip = self.nip + inner.nip

    def check_in(self, point):
        return None

class Circle(Shape):
    def __init__(self, nn, radius, center, domain):
        Shape.__init__(self, nn, nn, domain)
        self.radius = radius
        self.center = center
        
        self.nodes = [Node(center.x+radius*math.cos(i*2*math.pi/nn),center.y+radius*math.sin(i*2*math.pi/nn)) for i in range(nn)]
        Shape.gen_link(self)

    def check_in(self, point):
        r = point - self.center
        return math.sqrt(r.x**2 + r.y**2) < self.radius

class Square(Shape):
    def __init__(self, corner, lx, ly, nx, ny, domain):
        Shape.__init__(self, (nx+ny)*2, (nx+ny)*2, domain)
        dlx = lx/nx
        dly = ly/ny
        x0 = corner.x
        y0 = corner.y
        for i in range(nx):
            self.nodes.append(Point(x0+dlx*i,y0      ))
        x0 = corner.x+lx
        y0 = corner.y
        for i in range(ny):
            self.nodes.append(Point(x0      ,y0+dly*i))
        x0 = corner.x+lx
        y0 = corner.y+ly
        for i in range(nx):
            self.nodes.append(Point(x0-dlx*i,y0      ))
        x0 = corner.x
        y0 = corner.y+ly
        for i in range(ny):
            self.nodes.append(Point(x0      ,y0-dly*i))

        Shape.gen_link(self)

    def check_in(self, point):
        return False # todo

class Ellipse(Shape):
    def __init__(self, nn, r1, r2, center, theta, domain):
        Shape.__init__(self, nn, nn, domain)
        self.r1 = r1
        self.r2 = r2
        self.center = center
        self.theta = theta
        
        self.nodes = [Node(center.x+r1*math.cos(i*2*math.pi/n)*math.cos(theta)-r2*math.sin(i*2*math.pi/n)*math.sin(theta),center.y+r1*math.cos(i*2*math.pi/n)*math.sin(theta)+r2*math.sin(i*2*math.pi/n)*math.cos(theta)) for i in range(nn)]
        Shape.gen_link(self)

    def check_in(self, point):
        return False # todo

class Gear(Shape):
    def __init__(self, nn, r1, r2, center, n, theta0, domain):
        Shape.__init__(self, nn, nn, domain)
        __alpha = r2
        __beta = math.acos(r1/r2)
        theta0 = theta0*np.pi/180.0
        x = np.array([__alpha*cos(__beta*sin(n*0.5*(theta)))*cos((theta))  for theta in [2*np.pi*i/nn for i in range(nn)]])
        y = np.array([__alpha*cos(__beta*sin(n*0.5*(theta)))*sin((theta))  for theta in [2*np.pi*i/nn for i in range(nn)]])
        x_rotated = center.x +  x*np.cos(theta0) + y*np.sin(theta0)
        y_rotated = center.y -  x*np.sin(theta0) + y*np.cos(theta0)
        # x = np.array([center.x+__alpha*cos(__beta*sin(n*0.5*(theta)))*cos((theta))  for theta in [2*np.pi*i/nn for i in range(nn)]])
        # y = np.array([center.x+__alpha*cos(__beta*sin(n*0.5*(theta)))*sin((theta))  for theta in [2*np.pi*i/nn for i in range(nn)]])
        # x_rotated =  x*np.cos(theta0) + y*np.sin(theta0)
        # y_rotated = -x*np.sin(theta0) + y*np.cos(theta0)
        #self.nodes = [Node(center.x+__alpha*cos(__beta*sin(n*0.5*(theta)))*cos((theta)),center.x+__alpha*cos(__beta*sin(n*0.5*(theta)))*sin((theta))) for theta in [2*math.pi*i/nn for i in range(nn)]]
        self.nodes = [Node(x,y) for (x,y) in zip(x_rotated,y_rotated)]
        Shape.gen_link(self)

class Star(Shape):
    def __init__(self, r, center, nn, theta0, domain):
        # print(r,center,nn,theta0,domain)
        Shape.__init__(self, nn, nn, domain)
        theta0 = theta0*np.pi/180.0
        # x = np.array([__alpha*cos(__beta*sin(n*0.5*(theta)))*cos((theta))  for theta in [2*np.pi*i/nn for i in range(nn)]])
        x = np.array([r*(1.0+0.3*cos(5*theta))*cos(theta)/1.3  for theta in [2*np.pi*i/nn for i in range(nn)]])
        y = np.array([r*(1.0+0.3*cos(5*theta))*sin(theta)/1.3  for theta in [2*np.pi*i/nn for i in range(nn)]])
        x_rotated = center.x +  x*np.cos(theta0) + y*np.sin(theta0)
        y_rotated = center.y -  x*np.sin(theta0) + y*np.cos(theta0)
        self.nodes = [Node(x,y) for (x,y) in zip(x_rotated,y_rotated)]
        Shape.gen_link(self)

class Sine(Shape):
    def __init__(self, nn, x1, x2, y0, amp, domain):
        Shape.__init__(self, nn, nn-1, domain)
        ul = x2 - x1
        self.nodes = [Node(x1 + t, y0 + amp*np.sin(2.0*np.pi/ul*t)) for t in np.linspace(0.0, ul, nn)]
        Shape.gen_link(self, closed=False)
        
    def chech_in(self, point):
        return False

class Inner:
    def __init__(self, nip):
        self.nip = nip
        self.points = []

# Grid状に内点を設定する(corner:左下隅)
class Grid(Inner):
    def __init__(self, corner, lx, ly, nx, ny, mode):
        Inner.__init__(self, nx*ny)
        # 長方形[corner.x,corner.x+lx] \times [corner.y,corner.y+ly]のすべての辺に点が乗るように配置
        if mode==1:
            for i in range(nx):
                for j in range(ny):
                    self.points.append(Point(corner.x+lx/(nx-1)*i,corner.y+ly/(ny-1)*j))
        elif mode==2:
            for j in range(ny):
                for i in range(nx):
                    self.points.append(Point(corner.x+lx/(nx-1)*i,corner.y+ly/(ny-1)*j))
        # 長方形の右，上の辺の上には点が乗らない
        elif mode==3:
            for i in range(nx):
                for j in range(ny):
                    self.points.append(Point(corner.x+lx/nx*i,corner.y+ly/ny*j))
        elif mode==4:
            for j in range(ny):
                for i in range(nx):
                    self.points.append(Point(corner.x+lx/nx*i,corner.y+ly/ny*j))

class Line(Inner):
    def __init__(self, n, x, y, l, theta):
        Inner.__init__(self, n)
        dl = l/(n-1)
        theta = theta*np.pi/180
        self.points = [Point(x+dl*i*math.cos(theta),y+dl*i*math.sin(theta)) for i in range(n)]

class Square_in(Inner):
    def __init__(self, corner, lx, ly, nx, ny):
        Inner.__init__(self, (nx+ny)*2)
        dlx = lx/nx
        dly = ly/ny
        x0 = corner.x
        y0 = corner.y
        for i in range(nx):
            self.points.append(Point(x0+dlx*i,y0      ))
        x0 = corner.x+lx
        y0 = corner.y
        for i in range(ny):
            self.points.append(Point(x0      ,y0+dly*i))
        x0 = corner.x+lx
        y0 = corner.y+ly
        for i in range(nx):
            self.points.append(Point(x0-dlx*i,y0      ))
        x0 = corner.x
        y0 = corner.y+ly
        for i in range(ny):
            self.points.append(Point(x0      ,y0-dly*i))

class Circle_in(Shape):
    def __init__(self, nn, radius, center):
        Inner.__init__(self, nn)
        self.radius = radius
        self.center = center
        
        self.points = [Point(center.x+radius*math.cos(i*2*math.pi/n),center.y+radius*math.sin(i*2*math.pi/n)) for i in range(nn)]
                

canvas = Canvas()
# ファイルを読み込む
try:
    f = open("mesh.dat", "r")
except: # inputファイルがなければtemplateを出力して続行
    f = open("mesh.dat", "w")
    writer = csv.writer(f, lineterminator='\n',delimiter='\t')
    writer.writerow(["circle", 1000, 10.0, 0.0, 0.0, 2])
    writer.writerow(["naiten","grid", -20.0, -20.0, 40.0, 40.0, 80, 80])
    print("Wrote default input file 'mesh.dat'")
    f.close()
    f = open("mesh.dat", "r")

x1 = x2 = None
periodic = False
    
for line in f:
    #data = line[:-1].split("\t")
    data = line[:].split("\t")
    # 1列目がperiodicなら
    if data[0] == 'periodic':
        x1 = float(data[1])
        x2 = float(data[2])
        periodic = True
        
    # 1列目がcircleなら
    if data[0] == 'circle':
        n = int(data[1])
        radius = float(data[2])
        center = Point(float(data[3]),float(data[4]))
        domain = int(data[5])
        
        canvas.add(Circle(n,radius,center,domain))
    # 1列目がsquareなら
    if data[0] == 'square':
            corner = Point(float(data[1]),float(data[2]))
            lx = float(data[3])
            ly = float(data[4])
            nx = int(data[5])
            ny = int(data[6])
            domain = int(data[7])
            canvas.add(Square(corner,lx,ly,nx,ny,domain))
    # 1列目がellipseなら
    if data[0] == 'ellipse':
        n = int(data[1])
        r1 = float(data[2])
        r2 = float(data[3])
        center = Point(float(data[4]),float(data[5]))
        theta = float(data[6])*math.pi/180
        domain = int(data[7])
        canvas.add(Ellipse(n,r1,r2,center,theta,domain))
    # 1列目がgearなら
    if data[0] == 'gear':
        n = int(data[1])
        r1 = float(data[2])
        r2 = float(data[3])
        center = Point(float(data[4]),float(data[5]))
        n2 = int(data[6])
        theta0 = float(data[7])
        domain = int(data[8])
        canvas.add(Gear(n,r1,r2,center,n2,theta0,domain))
    # 1列目がstarなら
    if data[0] == 'star':
        n = int(data[1])
        r = float(data[2])
        center = Point(float(data[3]),float(data[4]))
        theta0 = float(data[5])
        domain = int(data[6])
        canvas.add(Star(r,center,n,theta0,domain))

    # 1列目がsineなら
    if data[0] == 'sine':
        n = int(data[1])
        y0 = float(data[2])
        amp = float(data[3])
        domain = int(data[4])
        canvas.add(Sine(n, x1, x2, y0, amp, domain))

    # 1列目がnaitenなら
    if data[0] == 'naiten':
        # 2列目がgridなら
        if data[1] == 'grid':
            corner = Point(float(data[2]),float(data[3]))
            lx = float(data[4])
            ly = float(data[5])
            nx = int(data[6])
            ny = int(data[7])
            mode = int(data[8])
            canvas.add_inner(Grid(corner,lx,ly,nx,ny,mode))
        if data[1] == 'line':
            n = int(data[2])
            x = float(data[3])
            y = float(data[4])
            l = float(data[5])
            theta = float(data[6])
            canvas.add_inner(Line(n,x,y,l,theta))
        if data[1] == 'square':
            corner = Point(float(data[2]),float(data[3]))
            lx = float(data[4])
            ly = float(data[5])
            nx = int(data[6])
            ny = int(data[7])
            canvas.add_inner(Square_in(corner,lx,ly,nx,ny))
        if data[1] == 'circle':
            n = int(data[2])
            radius = float(data[3])
            center = Point(float(data[4]),float(data[5]))
            canvas.add_inner(Circle_in(n,radius,center))

# ファイルに書き出す
if len(canvas.shapes) != 0:
    if periodic: filename = 'elm_init.el2p'
    else: filename = 'elm_init.el2'
    with open(filename, 'w') as f:
        writer = csv.writer(f, lineterminator='\n',delimiter='\t')
        if periodic:
            writer.writerow([x1,x2])
        # Nodeを書き込む
        writer.writerow([canvas.nn])
        for shape in canvas.shapes:
            for node in shape.nodes:
                row = [node.n+1, node.x, node.y]
                writer.writerow(row)
        # Linkを書き込む
        writer.writerow([canvas.nl])        
        for shape in canvas.shapes:
            random.shuffle(shape.links)
            if(shape.domain == 1):
                for link in shape.links:
                    row = [link.n+1, link.node2.n+1, link.node1.n+1]
                    writer.writerow(row)
            if(shape.domain == 2):
                for link in shape.links:
                    row = [link.n+1, link.node1.n+1, link.node2.n+1]
                    writer.writerow(row)        
    f.close()

if len(canvas.inners) != 0:
    with open('naiten.dat', 'w') as f:
        writer = csv.writer(f, lineterminator='\n',delimiter='\t')
        # 内点を書き込む
        writer.writerow([canvas.nip])
        i = 1
        for inner in canvas.inners:
            for point in inner.points:
                # id = 1
                # for shape in canvas.shapes:
                #     if(shape.check_in(point)): id = 2
                row = [i, point.x, point.y]
                writer.writerow(row)
                i = i+1
    f.close()
