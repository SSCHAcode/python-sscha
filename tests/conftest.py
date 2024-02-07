# -*- coding: utf-8 -*-
# pylint: disable=redefined-outer-name,too-many-statements
"""Initialise a text database and profile for pytest."""
from collections.abc import Mapping
import io
import os
import pathlib
import shutil
import tempfile

import pytest

pytest_plugins = ['aiida.manage.tests.pytest_fixtures']  # pylint: disable=invalid-name


@pytest.fixture(scope='session')
def filepath_tests():
    """Return the absolute filepath of the `tests` folder.

    .. warning:: if this file moves with respect to the `tests` folder, the implementation should change.

    :return: absolute filepath of `tests` folder which is the basepath for all test resources.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def filepath_fixtures(filepath_tests):
    """Return the absolute filepath to the directory containing the file `fixtures`."""
    return os.path.join(filepath_tests, 'fixtures')


@pytest.fixture(scope='function')
def fixture_sandbox():
    """Return a `SandboxFolder`."""
    from aiida.common.folders import SandboxFolder
    with SandboxFolder() as folder:
        yield folder


@pytest.fixture
def fixture_localhost(aiida_localhost):
    """Return a localhost `Computer`."""
    localhost = aiida_localhost
    localhost.set_default_mpiprocs_per_machine(1)
    return localhost


@pytest.fixture
def fixture_code(fixture_localhost):
    """Return an ``InstalledCode`` instance configured to run calculations of given entry point on localhost."""

    def _fixture_code(entry_point_name):
        from aiida.common import exceptions
        from aiida.orm import InstalledCode, load_code

        label = f'test.{entry_point_name}'

        try:
            return load_code(label=label)
        except exceptions.NotExistent:
            return InstalledCode(
                label=label,
                computer=fixture_localhost,
                filepath_executable='/bin/true',
                default_calc_job_plugin=entry_point_name,
            )

    return _fixture_code


@pytest.fixture
def serialize_builder():
    """Serialize the given process builder into a dictionary with nodes turned into their value representation.

    :param builder: the process builder to serialize
    :return: dictionary
    """

    def serialize_data(data):
        # pylint: disable=too-many-return-statements
        from aiida.orm import AbstractCode, BaseType, Data, Dict, KpointsData, List, RemoteData, SinglefileData
        from aiida.plugins import DataFactory

        StructureData = DataFactory('core.structure')
        UpfData = DataFactory('pseudo.upf')

        if isinstance(data, dict):
            return {key: serialize_data(value) for key, value in data.items()}

        if isinstance(data, BaseType):
            return data.value

        if isinstance(data, AbstractCode):
            return data.full_label

        if isinstance(data, Dict):
            return data.get_dict()

        if isinstance(data, List):
            return data.get_list()

        if isinstance(data, StructureData):
            return data.get_formula()

        if isinstance(data, UpfData):
            return f'{data.element}<md5={data.md5}>'

        if isinstance(data, RemoteData):
            # For `RemoteData` we compute the hash of the repository. The value returned by `Node._get_hash` is not
            # useful since it includes the hash of the absolute filepath and the computer UUID which vary between tests
            return data.base.repository.hash()

        if isinstance(data, KpointsData):
            try:
                return data.get_kpoints()
            except AttributeError:
                return data.get_kpoints_mesh()

        if isinstance(data, SinglefileData):
            return data.get_content()

        if isinstance(data, Data):
            return data.base.caching._get_hash()  # pylint: disable=protected-access

        return data

    def _serialize_builder(builder):
        return serialize_data(builder._inputs(prune=True))  # pylint: disable=protected-access

    return _serialize_builder


@pytest.fixture(scope='session', autouse=True)
def sssp(aiida_profile, generate_upf_data):
    """Create an SSSP pseudo potential family from scratch."""
    from aiida.common.constants import elements
    from aiida.plugins import GroupFactory

    aiida_profile.clear_profile()

    SsspFamily = GroupFactory('pseudo.family.sssp')

    cutoffs = {}
    stringency = 'standard'

    with tempfile.TemporaryDirectory() as dirpath:
        for values in elements.values():

            element = values['symbol']

            actinides = ('Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr')

            if element in actinides:
                continue

            upf = generate_upf_data(element)
            dirpath = pathlib.Path(dirpath)
            filename = dirpath / f'{element}.upf'

            with open(filename, 'w+b') as handle:
                with upf.open(mode='rb') as source:
                    handle.write(source.read())
                    handle.flush()

            cutoffs[element] = {
                'cutoff_wfc': 30.0,
                'cutoff_rho': 240.0,
            }

        label = 'SSSP/1.3/PBEsol/efficiency'
        family = SsspFamily.create_from_folder(dirpath, label)

    family.set_cutoffs(cutoffs, stringency, unit='Ry')

    return family


@pytest.fixture
def generate_calc_job():
    """Fixture to construct a new `CalcJob` instance and call `prepare_for_submission` for testing `CalcJob` classes.

    The fixture will return the `CalcInfo` returned by `prepare_for_submission` and the temporary folder that was passed
    to it, into which the raw input files will have been written.
    """

    def _generate_calc_job(folder, entry_point_name, inputs=None):
        """Fixture to generate a mock `CalcInfo` for testing calculation jobs."""
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager
        from aiida.plugins import CalculationFactory

        manager = get_manager()
        runner = manager.get_runner()

        process_class = CalculationFactory(entry_point_name)
        process = instantiate_process(runner, process_class, **inputs)

        calc_info = process.prepare_for_submission(folder)

        return calc_info

    return _generate_calc_job


@pytest.fixture
def generate_calc_job_node(fixture_localhost):
    """Fixture to generate a mock `CalcJobNode` for testing parsers."""

    def flatten_inputs(inputs, prefix=''):
        """Flatten inputs recursively like :meth:`aiida.engine.processes.process::Process._flatten_inputs`."""
        flat_inputs = []
        for key, value in inputs.items():
            if isinstance(value, Mapping):
                flat_inputs.extend(flatten_inputs(value, prefix=prefix + key + '__'))
            else:
                flat_inputs.append((prefix + key, value))
        return flat_inputs

    def _generate_calc_job_node(
        entry_point_name='base', computer=None, test_name=None, inputs=None, attributes=None, retrieve_temporary=None
    ):
        """Fixture to generate a mock `CalcJobNode` for testing parsers.

        :param entry_point_name: entry point name of the calculation class
        :param computer: a `Computer` instance
        :param test_name: relative path of directory with test output files in the `fixtures/{entry_point_name}` folder.
        :param inputs: any optional nodes to add as input links to the corrent CalcJobNode
        :param attributes: any optional attributes to set on the node
        :param retrieve_temporary: optional tuple of an absolute filepath of a temporary directory and a list of
            filenames that should be written to this directory, which will serve as the `retrieved_temporary_folder`.
            For now this only works with top-level files and does not support files nested in directories.
        :return: `CalcJobNode` instance with an attached `FolderData` as the `retrieved` node.
        """
        from aiida import orm
        from aiida.common import LinkType
        from aiida.plugins.entry_point import format_entry_point_string

        if computer is None:
            computer = fixture_localhost

        filepath_folder = None

        if test_name is not None:
            basepath = os.path.dirname(os.path.abspath(__file__))
            filename = os.path.join(entry_point_name[len('quantumespresso.'):], test_name)
            filepath_folder = os.path.join(basepath, 'parsers', 'fixtures', filename)
            filepath_input = os.path.join(filepath_folder, 'aiida.in')

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        node = orm.CalcJobNode(computer=computer, process_type=entry_point)
        node.base.attributes.set('input_filename', 'aiida.in')
        node.base.attributes.set('output_filename', 'aiida.out')
        node.base.attributes.set('error_filename', 'aiida.err')
        node.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
        node.set_option('max_wallclock_seconds', 1800)

        if attributes:
            node.base.attributes.set_many(attributes)

        if filepath_folder:
            from qe_tools.exceptions import ParsingError

            from aiida_quantumespresso.tools.pwinputparser import PwInputFile
            try:
                with open(filepath_input, 'r', encoding='utf-8') as input_file:
                    parsed_input = PwInputFile(input_file.read())
            except (ParsingError, FileNotFoundError):
                pass
            else:
                inputs['structure'] = parsed_input.get_structuredata()
                inputs['parameters'] = orm.Dict(parsed_input.namelists)

        if inputs:
            metadata = inputs.pop('metadata', {})
            options = metadata.get('options', {})

            for name, option in options.items():
                node.set_option(name, option)

            for link_label, input_node in flatten_inputs(inputs):
                input_node.store()
                node.base.links.add_incoming(input_node, link_type=LinkType.INPUT_CALC, link_label=link_label)

        node.store()

        if retrieve_temporary:
            dirpath, filenames = retrieve_temporary
            for filename in filenames:
                try:
                    shutil.copy(os.path.join(filepath_folder, filename), os.path.join(dirpath, filename))
                except FileNotFoundError:
                    pass  # To test the absence of files in the retrieve_temporary folder

        if filepath_folder:
            retrieved = orm.FolderData()
            retrieved.base.repository.put_object_from_tree(filepath_folder)

            # Remove files that are supposed to be only present in the retrieved temporary folder
            if retrieve_temporary:
                for filename in filenames:
                    try:
                        retrieved.base.repository.delete_object(filename)
                    except OSError:
                        pass  # To test the absence of files in the retrieve_temporary folder

            retrieved.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='retrieved')
            retrieved.store()

            remote_folder = orm.RemoteData(computer=computer, remote_path='/tmp')
            remote_folder.base.links.add_incoming(node, link_type=LinkType.CREATE, link_label='remote_folder')
            remote_folder.store()

        return node

    return _generate_calc_job_node


@pytest.fixture(scope='session')
def generate_upf_data():
    """Return a `UpfData` instance for the given element a file for which should exist in `tests/fixtures/pseudos`."""

    def _generate_upf_data(element):
        """Return `UpfData` node."""
        from aiida_pseudo.data.pseudo import UpfData
        content = f'<UPF version="2.0.1"><PP_HEADER\nelement="{element}"\nz_valence="4.0"\n/></UPF>\n'
        stream = io.BytesIO(content.encode('utf-8'))
        return UpfData(stream, filename=f'{element}.upf')

    return _generate_upf_data


@pytest.fixture(scope='session')
def generate_xy_data():
    """Return an ``XyData`` instance."""

    def _generate_xy_data():
        """Return an ``XyData`` node."""
        from aiida.orm import XyData
        import numpy as np

        xvals = [1, 2, 3]
        yvals = [10, 20, 30]
        xlabel = 'X'
        ylabel = 'Y'
        xunits = 'n/a'
        yunits = 'n/a'

        xy_node = XyData()
        xy_node.set_x(np.array(xvals), xlabel, xunits)
        xy_node.set_y([np.array(yvals)], [ylabel], [yunits])

        return xy_node

    return _generate_xy_data


@pytest.fixture
def generate_structure():
    """Return a ``StructureData`` representing either bulk silicon or a water molecule."""

    def _generate_structure(structure_id='silicon'):
        """Return a ``StructureData`` representing bulk silicon or a snapshot of a single water molecule dynamics.

        :param structure_id: identifies the ``StructureData`` you want to generate. Either 'silicon' or 'water'.
        """
        from aiida.orm import StructureData

        if structure_id.startswith('silicon'):
            name1 = 'Si0' if structure_id.endswith('kinds') else 'Si'
            name2 = 'Si1' if structure_id.endswith('kinds') else 'Si'
            param = 5.43
            cell = [[param / 2., param / 2., 0], [param / 2., 0, param / 2.], [0, param / 2., param / 2.]]
            structure = StructureData(cell=cell)
            structure.append_atom(position=(0., 0., 0.), symbols='Si', name=name1)
            structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='Si', name=name2)
        elif structure_id == 'water':
            structure = StructureData(cell=[[5.29177209, 0., 0.], [0., 5.29177209, 0.], [0., 0., 5.29177209]])
            structure.append_atom(position=[12.73464656, 16.7741411, 24.35076238], symbols='H', name='H')
            structure.append_atom(position=[-29.3865565, 9.51707929, -4.02515904], symbols='H', name='H')
            structure.append_atom(position=[1.04074437, -1.64320127, -1.27035021], symbols='O', name='O')
        elif structure_id == 'uranium':
            param = 5.43
            cell = [[param / 2., param / 2., 0], [param / 2., 0, param / 2.], [0, param / 2., param / 2.]]
            structure = StructureData(cell=cell)
            structure.append_atom(position=(0., 0., 0.), symbols='U', name='U')
            structure.append_atom(position=(param / 4., param / 4., param / 4.), symbols='U', name='U')
        elif structure_id == '2D-xy-arsenic':
            cell = [[3.61, 0, 0], [-1.80, 3.13, 0], [0, 0, 21.3]]
            structure = StructureData(cell=cell, pbc=(True, True, False))
            structure.append_atom(position=(1.804, 1.042, 11.352), symbols='As', name='As')
            structure.append_atom(position=(0, 2.083, 9.960), symbols='As', name='As')
        elif structure_id == '1D-x-carbon':
            cell = [[4.2, 0, 0], [0, 20, 0], [0, 0, 20]]
            structure = StructureData(cell=cell, pbc=(True, False, False))
            structure.append_atom(position=(0, 0, 0), symbols='C', name='C')
        elif structure_id == '1D-y-carbon':
            cell = [[20, 0, 0], [0, 4.2, 0], [0, 0, 20]]
            structure = StructureData(cell=cell, pbc=(False, True, False))
            structure.append_atom(position=(0, 0, 0), symbols='C', name='C')
        elif structure_id == '1D-z-carbon':
            cell = [[20, 0, 0], [0, 20, 0], [0, 0, 4.2]]
            structure = StructureData(cell=cell, pbc=(False, False, True))
            structure.append_atom(position=(0, 0, 0), symbols='C', name='C')
        else:
            raise KeyError(f'Unknown structure_id="{structure_id}"')
        return structure

    return _generate_structure


@pytest.fixture
def generate_kpoints_mesh():
    """Return a `KpointsData` node."""

    def _generate_kpoints_mesh(npoints):
        """Return a `KpointsData` with a mesh of npoints in each direction."""
        from aiida.orm import KpointsData

        kpoints = KpointsData()
        kpoints.set_kpoints_mesh([npoints] * 3)

        return kpoints

    return _generate_kpoints_mesh


@pytest.fixture(scope='session')
def generate_parser():
    """Fixture to load a parser class for testing parsers."""

    def _generate_parser(entry_point_name):
        """Fixture to load a parser class for testing parsers.

        :param entry_point_name: entry point name of the parser class
        :return: the `Parser` sub class
        """
        from aiida.plugins import ParserFactory
        return ParserFactory(entry_point_name)

    return _generate_parser


@pytest.fixture
def generate_remote_data():
    """Return a `RemoteData` node."""

    def _generate_remote_data(computer, remote_path, entry_point_name=None):
        """Return a `KpointsData` with a mesh of npoints in each direction."""
        from aiida.common.links import LinkType
        from aiida.orm import CalcJobNode, RemoteData
        from aiida.plugins.entry_point import format_entry_point_string

        entry_point = format_entry_point_string('aiida.calculations', entry_point_name)

        remote = RemoteData(remote_path=remote_path)
        remote.computer = computer

        if entry_point_name is not None:
            creator = CalcJobNode(computer=computer, process_type=entry_point)
            creator.set_option('resources', {'num_machines': 1, 'num_mpiprocs_per_machine': 1})
            remote.base.links.add_incoming(creator, link_type=LinkType.CREATE, link_label='remote_folder')
            creator.store()

        return remote

    return _generate_remote_data


@pytest.fixture
def generate_bands_data():
    """Return a `BandsData` node."""

    def _generate_bands_data():
        """Return a `BandsData` instance with some basic `kpoints` and `bands` arrays."""
        from aiida.plugins import DataFactory
        import numpy
        BandsData = DataFactory('core.array.bands')  #pylint: disable=invalid-name
        bands_data = BandsData()

        bands_data.set_kpoints(numpy.array([[0., 0., 0.], [0.625, 0.25, 0.625]]))

        bands_data.set_bands(
            numpy.array([[-5.64024889, 6.66929678, 6.66929678, 6.66929678, 8.91047649],
                         [-1.71354964, -0.74425095, 1.82242466, 3.98697455, 7.37979746]]),
            units='eV'
        )

        return bands_data

    return _generate_bands_data


@pytest.fixture
def generate_workchain():
    """Generate an instance of a `WorkChain`."""

    def _generate_workchain(entry_point, inputs):
        """Generate an instance of a `WorkChain` with the given entry point and inputs.

        :param entry_point: entry point name of the work chain subclass.
        :param inputs: inputs to be passed to process construction.
        :return: a `WorkChain` instance.
        """
        from aiida.engine.utils import instantiate_process
        from aiida.manage.manager import get_manager
        from aiida.plugins import WorkflowFactory

        process_class = WorkflowFactory(entry_point)
        runner = get_manager().get_runner()
        process = instantiate_process(runner, process_class, **inputs)

        return process

    return _generate_workchain


@pytest.fixture
def generate_force_constants_data(filepath_tests):
    """Generate a ``ForceConstantsData`` node."""
    from aiida_quantumespresso.data.force_constants import ForceConstantsData
    filepath = os.path.join(filepath_tests, 'calculations', 'fixtures', 'matdyn', 'default', 'force_constants.dat')
    return ForceConstantsData(filepath)


@pytest.fixture
def generate_inputs_pw(fixture_code, generate_structure, generate_kpoints_mesh, generate_upf_data):
    """Generate default inputs for a `PwCalculation."""

    def _generate_inputs_pw():
        """Generate default inputs for a `PwCalculation."""
        from aiida.orm import Dict

        from aiida_quantumespresso.utils.resources import get_default_options

        parameters = Dict({
            'CONTROL': {
                'calculation': 'scf'
            },
            'SYSTEM': {
                'ecutrho': 240.0,
                'ecutwfc': 30.0
            },
            'ELECTRONS': {
                'electron_maxstep': 60,
            }
        })
        structure = generate_structure()
        inputs = {
            'code': fixture_code('quantumespresso.pw'),
            'structure': generate_structure(),
            'kpoints': generate_kpoints_mesh(2),
            'parameters': parameters,
            'pseudos': {kind: generate_upf_data(kind) for kind in structure.get_kind_names()},
            'metadata': {
                'options': get_default_options()
            }
        }
        return inputs

    return _generate_inputs_pw


@pytest.fixture
def generate_workchain_pw_node():
    """Generate an instance of `WorkflowNode`."""
    from plumpy.process_states import ProcessState

    def _generate_workchain_pw_node(
        energy=-1.0,
        forces=None, 
        stress=None,
        label=None,
        exit_status=0, 
        process_state=ProcessState.FINISHED, 
    ):
        import numpy as np
        from aiida.common import LinkType
        from aiida.orm import TrajectoryData, Dict, WorkflowNode

        node = WorkflowNode()

        node.store()
        node.set_exit_status(exit_status)
        node.set_process_state(process_state)

        parameters = Dict({'energy': energy}).store()
        parameters.base.links.add_incoming(node, link_type=LinkType.RETURN, link_label='output_parameters')

        trajectory = TrajectoryData()
        if forces is None:
            trajectory.set_array('forces', np.zeros((1,3,2)))
        else:
            trajectory.set_array('forces', np.array(forces))
        if stress is None:
            trajectory.set_array('stress', np.zeros((1,3,3)))
        else:
            trajectory.set_array('stress', np.array(stress))
        
        stepids = np.array([1])
        times = stepids * 0.0
        cells = np.array([[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]])
        positions = np.array([[[0., 0., 0.], [0., 0., 0.]]])
        
        trajectory.set_trajectory(stepids=stepids, cells=cells, symbols=['Mg', 'O'], positions=positions, times=times)
        trajectory.store()
        trajectory.base.links.add_incoming(node, link_type=LinkType.RETURN, link_label='output_trajectory')
        
        if label:
            node.label = label

        return node

    return _generate_workchain_pw_node