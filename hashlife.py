import weakref
import math

class Universe:
    def round(self):
        """Compute (in place) the next generation of the universe"""
        raise NotImplementedError

    def get(self, i, j):
        """Returns the state of the cell at coordinates (ij[0], ij[1])"""
        raise NotImplementedError

    def rounds(self, n):
        """Compute (in place) the n-th next generation of the universe"""
        for _i in range(n):
            self.round()


class NaiveUniverse(Universe):
    def __init__(self, n, m, cells):
        self.n, self.m, self.cells = n, m, cells

    def round(self):
        new_cells = [[x for x in y] for y in self.cells]
        for i in range(self.n):
            for j in range(self.m):
                c = 0
                # check every cell around (i, j)
                if i > 0 and j > 0: c += 1 if self.cells[i - 1][j - 1] else 0
                if i > 0: c += 1 if self.cells[i - 1][j] else 0
                if j > 0: c += 1 if self.cells[i][j - 1] else 0
                if i > 0 and j + 1 < self.m: c += 1 if self.cells[i - 1][j + 1] else 0
                if j + 1 < self.m: c += 1 if self.cells[i][j + 1] else 0
                if i + 1 < self.n: c += 1 if self.cells[i + 1][j] else 0
                if i + 1 < self.n and j + 1 < self.m: c += 1 if self.cells[i + 1][j + 1] else 0
                if i + 1 < self.n and j > 0: c += 1 if self.cells[i + 1][j - 1] else 0
                new_cells[i][j] = c == 3 or c == 2 and self.cells[i][j]
        self.cells = new_cells

    def get(self, i, j):
        return self.cells[i][j] if 0 <= i < self.n and 0 <= j < self.m else False


class AbstractNode:
    BIG = True

    def __hash__(self):
        # lazy initialization
        res = getattr(self, "_hash", None)
        if res is None:
            self._hash = res = hash((
                self.population,
                self.level,
                self.nw,
                self.ne,
                self.sw,
                self.se,
            ))
        return res

    def __eq__(self, other):
        if self is other:
            return True
        if not isinstance(other, AbstractNode):
            return False
        return \
            self.level == other.level and \
            self.population == other.population and \
            self.nw is other.nw and \
            self.ne is other.ne and \
            self.sw is other.sw and \
            self.se is other.se

    @property
    def level(self):
        """Level of this node"""
        raise NotImplementedError()

    @property
    def population(self):
        """Total population of the area"""
        raise NotImplementedError()

    nw = property(lambda self: None)
    ne = property(lambda self: None)
    sw = property(lambda self: None)
    se = property(lambda self: None)

    @staticmethod
    def zero(k):
        res = AbstractNode.cell(0)
        for _ in range(k):
            res = AbstractNode.node(res, res, res, res)
        return res

    def extend(self):
        if isinstance(self, CellNode):
            return AbstractNode.node(
                AbstractNode.cell(0),
                self,
                AbstractNode.cell(0),
                AbstractNode.cell(0))
        zero = AbstractNode.zero(self.level - 1)
        return AbstractNode.node(
            AbstractNode.node(zero, zero, zero, self.nw),
            AbstractNode.node(zero, zero, self.ne, zero),
            AbstractNode.node(zero, self.sw, zero, zero),
            AbstractNode.node(self.se, zero, zero, zero)
        )

    _cache = None

    def forward(self, l=None):
        if self._cache is None: self._cache = {}
        l = self.level - 2 if l is None else min(max(l, 0), self.level - 2)
        if l in self._cache:
            return self._cache[l]
        if self.level < 2:
            s = None
        elif self.population < 3:
            s = AbstractNode.zero(self.level - 1)
        elif self.level == 2:
            w = self.nw.nw.population << 15 | self.nw.ne.population << 14 | \
                self.ne.nw.population << 13 | self.ne.ne.population << 12 | \
                self.nw.sw.population << 11 | self.nw.se.population << 10 | \
                self.ne.sw.population << 9 | self.ne.se.population << 8 | \
                self.sw.nw.population << 7 | self.sw.ne.population << 6 | \
                self.se.nw.population << 5 | self.se.ne.population << 4 | \
                self.sw.sw.population << 3 | self.sw.se.population << 2 | \
                self.se.sw.population << 1 | self.se.se.population
            s = Node.level2_bitmask(w)
        else:
            c1 = AbstractNode.node(self.nw.nw, self.nw.ne, self.nw.sw, self.nw.se).forward(l)
            c2 = AbstractNode.node(self.nw.ne, self.ne.nw, self.nw.se, self.ne.sw).forward(l)
            c3 = AbstractNode.node(self.ne.nw, self.ne.ne, self.ne.sw, self.ne.se).forward(l)
            c4 = AbstractNode.node(self.nw.sw, self.nw.se, self.sw.nw, self.sw.ne).forward(l)
            c5 = AbstractNode.node(self.nw.se, self.ne.sw, self.sw.ne, self.se.nw).forward(l)
            c6 = AbstractNode.node(self.ne.sw, self.ne.se, self.se.nw, self.se.ne).forward(l)
            c7 = AbstractNode.node(self.sw.nw, self.sw.ne, self.sw.sw, self.sw.se).forward(l)
            c8 = AbstractNode.node(self.sw.ne, self.se.nw, self.sw.se, self.se.sw).forward(l)
            c9 = AbstractNode.node(self.se.nw, self.se.ne, self.se.sw, self.se.se).forward(l)
            s = AbstractNode.node(
                AbstractNode.node(c1.se, c2.sw, c4.ne, c5.nw),
                AbstractNode.node(c2.se, c3.sw, c5.ne, c6.nw),
                AbstractNode.node(c4.se, c5.sw, c7.ne, c8.nw),
                AbstractNode.node(c5.se, c6.sw, c8.ne, c9.nw)
            ) if l < self.level - 2 else AbstractNode.node(
                AbstractNode.node(c1, c2, c4, c5).forward(l),
                AbstractNode.node(c2, c3, c5, c6).forward(l),
                AbstractNode.node(c4, c5, c7, c8).forward(l),
                AbstractNode.node(c5, c6, c8, c9).forward(l)
            )
        self._cache[l] = s
        return s

    @staticmethod
    def canon(node):
        return AbstractNode._node_cache.setdefault(node, node)

    @staticmethod
    def cell(alive):
        return AbstractNode.canon(CellNode(alive))

    @staticmethod
    def node(nw, ne, sw, se):
        return AbstractNode.canon(Node(nw, ne, sw, se))

    def get(self, i, j):
        raise NotImplementedError


AbstractNode._node_cache = weakref.WeakValueDictionary()


class CellNode(AbstractNode):
    def __init__(self, alive):
        super().__init__()
        self._alive = bool(alive)

    level = property(lambda self: 0)
    population = property(lambda self: int(self._alive))
    alive = property(lambda self: self._alive)

    def get(self, i, j):
        return i == j == 0 and self._alive


class Node(AbstractNode):
    @staticmethod
    def level2_bitmask(w):
        nw, ne, sw, se = [bin(x).count("1") for x in [w & 60128, w & 30064, w & 3758, w & 1879]]
        return AbstractNode.node(
            AbstractNode.cell(nw == 3 or nw == 2 and w & 1024),
            AbstractNode.cell(ne == 3 or ne == 2 and w & 512),
            AbstractNode.cell(sw == 3 or sw == 2 and w & 64),
            AbstractNode.cell(se == 3 or se == 2 and w & 32)
        )

    def __init__(self, nw, ne, sw, se):
        super().__init__()

        self._level = 1 + nw.level
        self._population = \
            nw.population + \
            ne.population + \
            sw.population + \
            se.population
        self._nw = nw
        self._ne = ne
        self._sw = sw
        self._se = se

    level = property(lambda self: self._level)
    population = property(lambda self: self._population)

    nw = property(lambda self: self._nw)
    ne = property(lambda self: self._ne)
    sw = property(lambda self: self._sw)
    se = property(lambda self: self._se)

    def get(self, i, j):
        lim = 2 ** self.level
        return [self.nw, self.ne, self.sw, self.se][(2 if i >= lim // 2 else 0) + (1 if j >= lim // 2 else 0)] \
            .get(i % (lim // 2), j % (lim // 2))


class HashLifeUniverse(Universe):
    def __init__(self, *args):
        if len(args) == 1:
            self._root = args[0]
        else:
            self._root = HashLifeUniverse.load(*args)

        self._generation = 0

    @staticmethod
    def load(n, m, cells):
        level = math.ceil(math.log(max(1, n, m), 2))

        mkcell = getattr(AbstractNode, 'cell', CellNode)
        mknode = getattr(AbstractNode, 'node', Node)

        def get(i, j):
            i, j = i + n // 2, j + m // 2
            return \
                i in range(n) and \
                j in range(m) and \
                cells[i][j]

        def create(i, j, level):
            if level == 0:
                return mkcell(get(i, j))

            noffset = 1 if level < 2 else 1 << (level - 2)
            poffset = 0 if level < 2 else 1 << (level - 2)

            nw = create(i - noffset, j + poffset, level - 1)
            sw = create(i - noffset, j - noffset, level - 1)
            ne = create(i + poffset, j + poffset, level - 1)
            se = create(i + poffset, j - noffset, level - 1)

            return mknode(nw=nw, ne=ne, sw=sw, se=se)

        return create(0, 0, level)

    def get(self, i, j):
        lim = 2 ** (self._root.level - 1)
        if i < -lim or i >= lim or j < -lim or j >= lim: return False
        node = self._root
        while True:
            if lim == 0:
                return node.population == 1
            node = [node.nw, node.ne, node.sw, node.se][(1 if i >= 0 else 0) + (0 if j >= 0 else 2)]
            lim //= 2
            i += lim * (int(i < 0) * 2 - 1)
            j += lim * (int(j < 0) * 2 - 1)

    def rounds(self, n):
        if n <= 0:
            raise ValueError("?")
        lvl = self.root.level
        k = 0
        while n > 0:
            self.extend(max(k + 2, self._root.level + 2))
            if n & 1:
                self._root = self._root.forward(k)
            k += 1
            n >>= 1
            center = AbstractNode.node(
                self._root.nw.se,
                self._root.ne.sw,
                self._root.sw.ne,
                self._root.se.nw
            )
            while center.population == self._root.population and self._root.level > 1:
                self._root = center
                center = AbstractNode.node(
                    self._root.nw.se,
                    self._root.ne.sw,
                    self._root.sw.ne,
                    self._root.se.nw
                )
        self._generation += n

    def round(self):
        return self.rounds(1)

    @property
    def root(self):
        return self._root

    @property
    def generation(self):
        return self._generation

    def extend(self, k):
        res = self._root
        # check the level and the peripheral band
        while res.level < max(k, 2) or \
                res.nw.nw.population > 0 or res.nw.ne.population > 0 or \
                res.ne.nw.population > 0 or res.ne.ne.population > 0 or \
                res.nw.sw.population > 0 or res.ne.se.population > 0 or \
                res.sw.nw.population > 0 or res.se.ne.population > 0 or \
                res.sw.sw.population > 0 or res.sw.se.population > 0 or \
                res.se.sw.population > 0 or res.se.se.population > 0:
            res = res.extend()
        self._root = res
