import json


class JsonStorage():
    def __init__(self, fullpath):
        self.fullpath = fullpath

    def get(self, name):
        with open(self.fullpath) as f:
            val = json.load(f)[name]
        return val
    
    def get_all(self):
        with open(self.fullpath) as f:
            return json.load(f)

    def save(self, name, val):
        all_new = self.get_all()
        all_new[name] = val
        with open(self.fullpath, 'w') as f:
            json.dump(all_new, f, indent=4)
