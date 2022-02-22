import cProfile
from dataclasses import dataclass
import math
import random

@dataclass(frozen=False)
class bhVector:
    x: float
    y: float
    z: float
    
    def length(self) -> float:
        return math.sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
    
    def normalize(self) -> None:
        length = self.length()
        self.x = self.x / length
        self.y = self.y / length
        self.z = self.z / length
        
    def multScalar(self, s: float) -> None:
        self.x *= s
        self.y *= s
        self.z *= s
        
    def add(v1, v2):
        return bhVector(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z)

@dataclass
class bhBody:
    x: float
    y: float
    z: float
    mass: float
    force: bhVector

class bhRegion:
    def __init__(self, x: float, y: float, z: float, w: float, h: float, d: float) -> None:
        # Center coordinates
        self.x = x
        self.y = y
        self.z = z
        self.w = w  # width
        self.h = h  # height
        self.d = d  # depth

    def contains(self, body: bhBody) -> bool:
        return (body.x <= self.x + .5 * self.w and
                body.x >= self.x - .5 * self.w and
                body.y <= self.y + .5 * self.h and
                body.y >= self.y - .5 * self.h and
                body.z <= self.z + .5 * self.d and
                body.z >= self.z - .5 * self.d)

    def __repr__(self) -> str:
        rep = f'x:{self.x}, y:{self.y}, z:{self.z}, w:{self.w}, h:{self.h}, d:{self.d}'
        return rep

class bhTree:
    def __init__(self, region: bhRegion, capacity: int, depth: int, theta: float) -> None:
        self.region = region
        self.capacity = capacity
        self.bodies = []
        self.external = True
        self.depth = depth
        self.totalMass = 0
        self.massCenter = bhBody(
            self.region.x, self.region.y, self.region.z, 0, bhVector(0, 0, 0))
        self.theta = theta # Parameter to eveluate whether we need to calculate with center of mass or continue in children node

    def subdivide(self) -> None:
        # +x towards the right / +y towards the bottom / +z towards the screen (user)
        fnw = bhRegion(self.region.x - self.region.w/4, self.region.y - self.region.h/4,
                       self.region.z + self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.frontnorthwest = bhTree(fnw, self.capacity, self.depth+1, self.theta)
        fne = bhRegion(self.region.x + self.region.w/4, self.region.y - self.region.h/4,
                       self.region.z + self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.frontnortheast = bhTree(fne, self.capacity, self.depth+1, self.theta)
        fsw = bhRegion(self.region.x - self.region.w/4, self.region.y + self.region.h/4,
                       self.region.z + self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.frontsouthwest = bhTree(fsw, self.capacity, self.depth+1, self.theta)
        fse = bhRegion(self.region.x + self.region.w/4, self.region.y + self.region.h/4,
                       self.region.z + self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.frontsoutheast = bhTree(fse, self.capacity, self.depth+1, self.theta)
        bnw = bhRegion(self.region.x - self.region.w/4, self.region.y - self.region.h/4,
                       self.region.z - self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.backnorthwest = bhTree(bnw, self.capacity, self.depth+1, self.theta)
        bne = bhRegion(self.region.x + self.region.w/4, self.region.y - self.region.h/4,
                       self.region.z - self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.backnortheast = bhTree(bne, self.capacity, self.depth+1, self.theta)
        bsw = bhRegion(self.region.x - self.region.w/4, self.region.y + self.region.h/4,
                       self.region.z - self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.backsouthwest = bhTree(bsw, self.capacity, self.depth+1, self.theta)
        bse = bhRegion(self.region.x + self.region.w/4, self.region.y + self.region.h/4,
                       self.region.z - self.region.d/4, self.region.w/2, self.region.h/2, self.region.d/2)
        self.backsoutheast = bhTree(bse, self.capacity, self.depth+1, self.theta)
        self.external = False
        for b in self.bodies:
            if (self.frontnorthwest.insert(b)):
                continue
            if (self.frontnortheast.insert(b)):
                continue
            if (self.frontsouthwest.insert(b)):
                continue
            if (self.frontsoutheast.insert(b)):
                continue
            if (self.backnorthwest.insert(b)):
                continue
            if (self.backnortheast.insert(b)):
                continue
            if (self.backsouthwest.insert(b)):
                continue
            if (self.backsoutheast.insert(b)):
                continue
        del self.bodies[:]

    def insert(self, body: bhBody) -> bool:
        if not (self.region.contains(body)):
            return False
        if (len(self.bodies) < self.capacity and self.external):
            self.totalMass += body.mass
            self.bodies.append(body)
            return True
        else:
            if (self.external):
                self.subdivide()
            if (self.frontnorthwest.insert(body)):
                return True
            if (self.frontnortheast.insert(body)):
                return True
            if (self.frontsouthwest.insert(body)):
                return True
            if (self.frontsoutheast.insert(body)):
                return True
            if (self.backnorthwest.insert(body)):
                return True
            if (self.backnortheast.insert(body)):
                return True
            if (self.backsouthwest.insert(body)):
                return True
            if (self.backsoutheast.insert(body)):
                return True

    def calculateMassCenter(self) -> bhBody:
        nullForce = bhVector(0, 0, 0)
        c_x, c_y, c_z, c_mass = 0, 0, 0, 0
        if not self.external:
            massBodies = []
            massBodies.append(self.frontnorthwest.calculateMassCenter())
            massBodies.append(self.frontnortheast.calculateMassCenter())
            massBodies.append(self.frontsouthwest.calculateMassCenter())
            massBodies.append(self.frontsoutheast.calculateMassCenter())
            massBodies.append(self.backnorthwest.calculateMassCenter())
            massBodies.append(self.backnortheast.calculateMassCenter())
            massBodies.append(self.backsouthwest.calculateMassCenter())
            massBodies.append(self.backsoutheast.calculateMassCenter())
            if not massBodies:
                return None
            for mB in massBodies:
                c_x += mB.x * mB.mass
                c_y += mB.y * mB.mass
                c_z += mB.z * mB.mass
                c_mass += mB.mass
            c_x /= c_mass
            c_y /= c_mass
            c_z /= c_mass
            self.massCenter = bhBody(c_x, c_y, c_z, c_mass, nullForce)
            return self.massCenter
        if self.bodies:
            for mB in self.bodies:
                c_x += mB.x * mB.mass
                c_y += mB.y * mB.mass
                c_z += mB.z * mB.mass
                c_mass += mB.mass
            c_x /= c_mass
            c_y /= c_mass
            c_z /= c_mass
            self.massCenter = bhBody(c_x, c_y, c_z, c_mass, nullForce)
            return self.massCenter
        return bhBody(self.region.x, self.region.y, self.region.z, 0, nullForce)

    def queryForceOnBody(self, body: bhBody) -> bool:
        sqDist = self.squaredDistanceBetweenBodies(body, self.massCenter)
        if (sqDist == 0):
            return True
        if (self.external):
            if (len(self.bodies) == 0): return True # Node is external and empty we return it has been processed
            self.calculateForceOnBody(body, self.massCenter, sqDist)
            return True
        if (self.region.w / math.sqrt(sqDist) < self.theta):
            self.calculateForceOnBody(body, self.massCenter, sqDist)
            return True
        self.frontnorthwest.queryForceOnBody(body)
        self.frontnortheast.queryForceOnBody(body)
        self.frontsouthwest.queryForceOnBody(body)
        self.frontsoutheast.queryForceOnBody(body)
        self.backnorthwest.queryForceOnBody(body)
        self.backnortheast.queryForceOnBody(body)
        self.backsouthwest.queryForceOnBody(body)
        self.backsoutheast.queryForceOnBody(body)
        return True

    def squaredDistanceBetweenBodies(self, b1: bhBody, b2: bhBody) -> float:
        return (b1.x - b2.x) ** 2 + (b1.y - b2.y) ** 2 + (b1.z - b2.z) ** 2
    
    def calculateForceOnBody(self, focusBody: bhBody, targetBody: bhBody, sqDist: float) -> None:
        # focusBody is the body on which the force will be exerted
        # targetBody is the body that exerts the force on focusBody
        # dist is the distance between the bodies, passed as arg because already calculated (to reduce calculations)
        gravityConst = 1
        force = gravityConst * (focusBody.mass * targetBody.mass) / sqDist # Newton's law of universal gravitation
        forceVect = bhVector(targetBody.x - focusBody.x, targetBody.y - focusBody.y, targetBody.z - focusBody.z)
        forceVect.normalize()
        forceVect.multScalar(force)
        focusBody.force =  bhVector.add(focusBody.force, forceVect)
        

    def __repr__(self) -> str:
        depthIndent = '\t'*self.depth
        rep = f'{self.depth}_{depthIndent}bhTree({self.region.__repr__()}, Bodies count: {len(self.bodies)}, Bodies pos: {self.bodies.__repr__()}, Mass center : {self.massCenter.__repr__()})'
        if not self.external:
            rep = f'{rep} \nSubTrees : \nFNW_{self.frontnorthwest.__repr__()}\nFNE_{self.frontnortheast.__repr__()}\nFSW_{self.frontsouthwest.__repr__()}\nFSE_{self.frontsoutheast.__repr__()}\nBNW_{self.backnorthwest.__repr__()}\nBNE_{self.backnortheast.__repr__()}\nBSW_{self.backsouthwest.__repr__()}\nBSE_{self.backsoutheast.__repr__()}\n '
        return rep

def main(globalBodies):
    tree = bhTree(bhRegion(0, 0, 0, 100, 100, 100), 2, 0, 0.5)
    bodies = []
    for b in globalBodies:
        bodies.append(b)
        tree.insert(b)
    tree.calculateMassCenter()
    for b in bodies:
        tree.queryForceOnBody(b)

# Test used to compare execution speed between Brute Force and Barnes-Hut algorithm
def testForce(globalBodies):
    for b in globalBodies:
        fVect = bhVector(0, 0, 0)
        for i in range(len(globalBodies)):
            if globalBodies[i] == b: continue
            force = 1 * (b.mass * globalBodies[i].mass) / ((b.x - globalBodies[i].x) ** 2 + (b.y - globalBodies[i].y) ** 2 + (b.z - globalBodies[i].z) ** 2)
            iVect = bhVector(globalBodies[i].x - b.x, globalBodies[i].y - b.y, globalBodies[i].z - b.z)
            iVect.normalize()
            iVect.multScalar(force)
            fVect = bhVector.add(fVect, iVect)


if __name__ == '__main__':
    bodies = []
    for i in range(2000):
        b = bhBody(random.uniform(-50, 50), random.uniform(-50, 50),
                   random.uniform(-50, 50), 1.0, bhVector(0, 0, 0))
        bodies.append(b)
    cProfile.runctx('main(bodies)', {'main': main, 'bodies': bodies}, {})
    cProfile.runctx('testForce(bodies)', {'testForce': testForce, 'bodies': bodies}, {})
