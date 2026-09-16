class Morphism:
    def __init__(self):
        self._domain = None
        self._codomain = None

    def domain(self):
        """Return the domain of the morphism"""
        return self._domain

    def codomain(self):
        """Return the codomain of the morphism"""
        return self._codomain

    def __call__(self, point):
        raise NotImplementedError("subclasses must implement morphism evaluation")
