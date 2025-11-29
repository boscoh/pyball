"""
Alternative approach: Wrapper classes that provide persistent object interface
while working with transient proxies under the hood.
"""

from pdbstruct import parse
from pdbstruct.vector3d import Vector3d

class AtomWrapper:
    """Wrapper around atom index that behaves like a persistent object"""
    def __init__(self, soup, atom_idx, rendered_soup):
        self._soup = soup
        self._atom_idx = atom_idx
        self._rendered_soup = rendered_soup
    
    @property
    def pos(self):
        """Get position from underlying proxy"""
        return self._soup.get_atom_proxy(self._atom_idx).pos
    
    @property
    def atom_type(self):
        return self._soup.get_atom_proxy(self._atom_idx).atom_type
    
    @property
    def element(self):
        return self._soup.get_atom_proxy(self._atom_idx).element
    
    @property
    def type(self):
        """Alias for atom_type for backwards compat"""
        return self.atom_type
    
    @property
    def objid(self):
        """Get objid from rendered soup storage"""
        return self._rendered_soup.atom_objids.get(self._atom_idx)
    
    @objid.setter
    def objid(self, value):
        self._rendered_soup.atom_objids[self._atom_idx] = value
    
    @property  
    def residue(self):
        """Get residue wrapper"""
        res_idx = self._rendered_soup.atom_residue_idx.get(self._atom_idx)
        if res_idx is not None:
            return ResidueWrapper(self._soup, res_idx, self._rendered_soup)
        return None
    
    @residue.setter
    def residue(self, res_wrapper):
        if isinstance(res_wrapper, ResidueWrapper):
            self._rendered_soup.atom_residue_idx[self._atom_idx] = res_wrapper._res_idx


class ResidueWrapper:
    """Wrapper around residue index that behaves like a persistent object"""
    def __init__(self, soup, res_idx, rendered_soup):
        self._soup = soup
        self._res_idx = res_idx
        self._rendered_soup = rendered_soup
    
    @property
    def res_type(self):
        return self._soup.get_residue_proxy(self._res_idx).res_type
    
    @property
    def res_num(self):
        return self._soup.get_residue_proxy(self._res_idx).res_num
    
    @property
    def chain(self):
        return self._soup.get_residue_proxy(self._res_idx).chain
    
    def get_atom_indices(self):
        return self._soup.get_residue_proxy(self._res_idx).get_atom_indices()
    
    def has_atom(self, atom_type):
        return self._rendered_soup.has_atom_in_residue_idx(self._res_idx, atom_type)
    
    def atom(self, atom_type):
        """Get atom wrapper by type"""
        atom_idx = self._rendered_soup.find_atom_in_residue_idx(self._res_idx, atom_type)
        if atom_idx is not None:
            return AtomWrapper(self._soup, atom_idx, self._rendered_soup)
        return None
    
    def atoms(self):
        """Iterator over atom wrappers"""
        for atom_idx in self.get_atom_indices():
            yield AtomWrapper(self._soup, atom_idx, self._rendered_soup)
    
    # Custom attributes stored in rendered_soup
    @property
    def ss(self):
        return self._rendered_soup.residue_ss.get(self._res_idx, '-')
    
    @ss.setter
    def ss(self, value):
        self._rendered_soup.residue_ss[self._res_idx] = value
    
    @property
    def color(self):
        return self._rendered_soup.residue_color.get(self._res_idx, [0.4, 1.0, 0.4])
    
    @color.setter
    def color(self, value):
        self._rendered_soup.residue_color[self._res_idx] = value
    
    @property
    def objid(self):
        return self._rendered_soup.residue_objids.get(self._res_idx)
    
    @objid.setter
    def objid(self, value):
        self._rendered_soup.residue_objids[self._res_idx] = value
    
    @property
    def i(self):
        return self._rendered_soup.residue_i.get(self._res_idx)
    
    @i.setter
    def i(self, value):
        self._rendered_soup.residue_i[self._res_idx] = value
    
    @property
    def hb_partners(self):
        return self._rendered_soup.residue_hb_partners.get(self._res_idx, [])
    
    @hb_partners.setter
    def hb_partners(self, value):
        self._rendered_soup.residue_hb_partners[self._res_idx] = value

# Test it
if __name__ == "__main__":
    from pathlib import Path
    pdb_file = str(Path(__file__).parent.parent / "1be9.pdb")
    
    soup = parse.load_soup(pdb_file)
    
    class MockRenderedSoup:
        def __init__(self):
            self.atom_objids = {}
            self.atom_residue_idx = {}
            self.residue_ss = {}
            self.residue_color = {}
            self.residue_objids = {}
            self.residue_i = {}
            self.residue_hb_partners = {}
        
        def has_atom_in_residue_idx(self, res_idx, atom_type):
            proxy = soup.get_residue_proxy(res_idx)
            for atom_idx in proxy.get_atom_indices():
                if soup.get_atom_proxy(atom_idx).atom_type == atom_type:
                    return True
            return False
        
        def find_atom_in_residue_idx(self, res_idx, atom_type):
            proxy = soup.get_residue_proxy(res_idx)
            for atom_idx in proxy.get_atom_indices():
                if soup.get_atom_proxy(atom_idx).atom_type == atom_type:
                    return atom_idx
            return None
    
    rs = MockRenderedSoup()
    
    # Test residue wrapper
    res = ResidueWrapper(soup, 0, rs)
    print(f"Residue: {res.res_type}{res.res_num} chain {res.chain}")
    
    # Test setting attributes
    res.ss = 'H'
    res.color = [1, 0, 0]
    print(f"SS: {res.ss}, Color: {res.color}")
    
    # Test atom access
    if res.has_atom('CA'):
        ca = res.atom('CA')
        print(f"CA pos: {ca.pos}")
        ca.objid = 42
        print(f"CA objid: {ca.objid}")
    
    print("\\n✅ Wrapper approach works!")

