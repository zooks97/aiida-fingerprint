import json
import flask
import flask_restful
from flask_restful.reqparse import Argument, RequestParser

import ase.io

from sys import path
path.insert(0, '/root/git/glosim2')
from libmatch.soap import get_soap, get_Soaps
from libmatch.utils import ase2qp, get_spkit, get_spkitMax

app = flask.Flask(__name__)
api = flask_restful.Api(app)
VERSION = 1

ARGUMENTS = {
    'atoms': Argument('atoms', type=str, required=True),
    'spkitMax': Argument('spkitMax', type=str, required=True),
    'spkit': Argument('spkit', type=str, default=None, required=False),
    'nocenters': Argument('nocenters', default=None, required=False),
    'centerweight': Argument('centerweight', type=float, default=1.0, required=False),
    'gaussian_width': Argument('gaussian_width', type=float, default=0.5, required=False),
    'cutoff': Argument('cutoff', type=float, default=3.5, required=False),
    'cutoff_transition_width': Argument('cutoff_transition_width', type=float, default=0.5, required=False),
    'nmax': Argument('nmax', type=int, default=8, required=False),
    'lmax': Argument('lmax', type=int, default=6, required=False),
    # 'is_fast_average': Argument('is_fast_average', type=bool, default=False, required=False),
}


class FakeFile(object):
    def __init__(self, contents):
        self.contents = contents

    def seek(self, arg):
        return

    def read(self):
        return self.contents


class soapv1(flask_restful.Resource):
    def get(self):
        parser = flask_restful.reqparse.RequestParser()
        argument_names = ['atoms', 'spkitMax', 'nocenters', 'gaussian_width', 'spkit',
                          'cutoff', 'cutoff_transition_width', 'nmax', 'lmax']
        for argument_name in argument_names:
            parser.add_argument(ARGUMENTS[argument_name])
        args = parser.parse_args(strict=True)
        args['atoms'] = ase.io.read(FakeFile(args['atoms']), format='json')
	args['atoms'] = ase2qp(args['atoms'])
        args['spkitMax'] = json.loads(args['spkitMax'])
        if not args['spkit']:
            args['spkit'] = get_spkit(args['atoms'])
	# args['atoms'] = [args['atoms']]
        soaps = get_soap(**args)
        # soaps = {key: value.tolist() for key, value in soaps.iteritems()}
	soaps = soaps.tolist()

        return flask.jsonify(soaps)


api.add_resource(soapv1, '/v{}/soap/'.format(VERSION))

if __name__ == '__main__':
    app.run(debug=True, port=8899, host='0.0.0.0')
