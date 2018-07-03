import flask
import flask_restful
from flask_restful.reqparse import Argument, RequestParser

import ase.io

path.insert(0, '/root/git/glosim2')
from libmatch.soap import get_soap
from libmatch.utils import ase2qp, get_spkit, get_spkitMax

app = flask.Flask(__name__)
api = flask_restful.Api(app)
VERSION = 1

ARGUMENTS = {
    'atoms': Argument('atoms', type=str, required=True, help='JSON string of quippy atoms object'),
    'spkit': Argument('spkit', type=dict, required=True, help='Dictionary of atomic numbers and number of occurences for the structure'),
    'spkitMax': Argument('spkitMax', type=dict, required=True, help='Dictionary of atomic numbers and maximum number of occurences in all structures'),
    'nocenters': Argument('nocenters', default=None),
    'centerweight': Argument('centerweight', type=float, default=1.0),
    'gaussian_width': Argument('gaussian_width', type=float, default=0.5),
    'cutoff': Argument('cutoff', type=float, default=3.5),
    'cutoff_transition_width': Argument('cutoff_transition_width', type=float, default=0.5),
    'nmax': Argument('nmax', type=int, default=8),
    'lmax': Argument('lmax', type=int, default=6),
    'is_fast_average': Argument('is_fast_average', type=bool, default=False),
}


class FakeFile(object):
    def __init__(self, contents):
        self.contents = contents

    def seek(self, arg):
        return

    def read(self)
        return self.contents


class soapv1(flask_restful.Resource):
    def get(self):
        parser = flask_restful.reqparse.RequestParser()
        argument_names = ['atoms', 'spkit', 'spkitMax', 'nocenters', 'gaussian_width',
                          'cutoff', 'cutoff_transition_width', 'nmax', 'lmax', 'is_fast_average']
        for argument_name in argument_names:
            parser.add_argument(ARGUMENTS[argument_name])
        args = parser.parse_args(strict=True)
        args['atoms'] = ase.io.read(FakeFile(args['atoms']), format='json')
        soaps = get_soap(**args)

        return flask.jsonify(dict(soaps))


api.add_resource(soapv1, 'v{}/soap/'.format(VERSION))
