"""Defines a round of the campaign.

`NoAideRound` representes a minimal viable round of the campaign, where all unlabeled
variants are selected for experiment. 

To modify the generation or selection strategy for a round, create a new mixin
class that overloads the _setup_library_generation or _setup_library_selection, and define
a new parameter class for the round. See GenerationRoundMixin and SelectionRoundMixin.

Eg. to define a new round within the AIDE paradigm from scratch, the recommended approach is:

```
from aide.campaign.round import BaseRound, RoundParams, GenerationRoundMixin, SelectionRoundMixin

class MyRoundParams(RoundParams):
    my_param: str = "my_default_value"

class MySelectionRoundMixin(SelectionRoundMixin):
    def _setup_library_selection(self) -> AquisitionFunction:
        # define a new selection strategy using self.params
        # must return a class that inherits from AquisitionFunction
        return my_selection_strategy

    def _any_other_method(self):
        # define helper functions for the above

class MyRound(GenerationRoundMixin, MySelectionRoundMixin, BaseRound):
    # here we left generation strategy as default, but could also overload similar to
    # MySelectionRoundMixin
    param_class = MyRoundParams
"""

from __future__ import annotations

from datetime import datetime

from aide.base import Library, EmptyLibraryException
from aide.base import LibraryGenerator, AquisitionFunction
from aide.campaign.database import CampaignDatabase
from dataclasses import dataclass

class RoundStatusError(Exception):
    """Raised when attempting to perform an action on a round with an invalid status."""
    pass

class RoundStatusWarning(Warning):
    """Raised when attempting to perform an action on a round with an invalid status."""
    pass

@dataclass
class RoundParams:
    _placeholder: str = "placeholder"

class BaseRound:
    """Abstract parent class for a round of the campaign.
    
    A round is a single iteration of the campaign, where a library is generated
    and then selected from.

    Subclass and define both _setup_library_generation and _setup_library_selection
    to attain functionaility. For a minimal viable round, use the Round class, which
    by default uses the AllUnlabeledGenerator and AllAquisition selection strategies
    eg. selects all variants from the database and does not downselect at all.

    Params
    ------
    database: CampaignDatabase
        Database to store and retrieve campaign data.
    params: RoundParams
        Parameters for the round. Must be of type RoundParams. For subclasses, create
        a new subclass of RoundParams to define the parameters for the round.

    Methods
    -------
    add_to_db
        Add the round to the database, and set the round index.
    commit
        Commit the round to the database. The round becomes immutable after this, unless reset.
    reset
        Reset the round to its initial state.
    library_generation
        Generate a library of variants, store in database, and setup for selection.
    library_selection
        Select a library from the putative library, and return the library for experiment.
    from_db
        Load a round from the database, given the database and round index.
    
    Attributes
    ----------
    db: CampaignDatabase
        Database to store and retrieve campaign data.
    round_idx: int
        Index of the round in the database.
    params: RoundParams
        Parameters for the round.
    putative_library: Library
        Library of variants generated for selection. Pupulated after library_generation.
    library_for_exp: Library
        Library of variants selected for experiment. Pupulated after library_selection.
    library_generator: LibraryGenerator
        Library generator for the round. Pupulated after _setup_library_generation.
    aquisition_function: AquisitionFunction
        Aquisition function for the round. Pupulated after _setup_library_selection.
    status: str
        Status of the round. One of "unknown", "not ready", "ready", "generated", "selected", "labeled", "complete".
    """
    param_class = None

    def __init__(
        self,
        database: CampaignDatabase,
        params: RoundParams=RoundParams(),
    ) -> BaseRound:
        self._database = database

        if not isinstance(params, self.param_class):
            raise ValueError(f"params must be of type {self.param_class.__name__}")
        self.params = params
        self._round_idx = None

        self.putative_library = None
        self.library_for_exp = None

        self.library_generator = self._setup_library_generation()
        if not isinstance(self.library_generator, LibraryGenerator):
            raise ValueError(f"library_generator must be of type LibraryGenerator")
        self.aquisition_function = self._setup_library_selection()
        if not isinstance(self.aquisition_function, AquisitionFunction):
            raise ValueError(f"aquisition_function must be of type AquisitionFunction")

    @property
    def db(self) -> CampaignDatabase:
        """Database to store and retrieve campaign data."""
        return self._database

    @property
    def round_idx(self) -> int:
        """Index of the round in the database."""
        return self._round_idx
    
    @property
    def status(self) -> str:
        """Status of the round. 
        
        One of "unknown", "not ready", "ready", "generated", "selected", "labeled", "complete".
        
        - unknown: round has not been added to database
        - not ready: round has been added to database, but not ready for library generation, eg. the previous round has not been committed
        - ready: round is ready for library generation
        - generated: putative library has been generated
        - selected: library for experiment has been selected
        - labeled: labels have been set for the library for experiment
        - complete: round has been locked to database, next round can be generated
        """
        if self._round_idx is None:
            return "unknown"
        else:
            return self.db.get_round_status(self.round_idx)
        
    @status.setter
    def status(self, status: str):
        if self._round_idx is None:
            raise RoundStatusError("Cannot set status for round that has not been committed to database.")
        self.db.set_round_status(self.round_idx, status)

    @property
    def notes(self) -> str:
        """Any experimental notes for the round."""
        if self._round_idx is None:
            return None
        else:
            return self.db.get_round_notes(self.round_idx)
        
    @notes.setter
    def notes(self, notes: str):
        if self._round_idx is None:
            raise RoundStatusError("Cannot set notes for round that has not been committed to database.")
        self.db.set_round_notes(self.round_idx, notes)

    @property
    def start_time(self) -> str:
        """Start time of the round, in format %Y-%m-%d %H:%M:%S."""
        if self._round_idx is None:
            return None
        else:
            return self.db.get_round_start_time(self.round_idx)
        
    @start_time.setter
    def start_time(self, start_time: str):
        if self._round_idx is None:
            raise RoundStatusError("Cannot set start time for round that has not been committed to database.")
        self.db.set_round_start_time(self.round_idx, start_time)

    @property
    def end_time(self) -> str:
        """End time of the round, in format %Y-%m-%d %H:%M:%S."""
        if self._round_idx is None:
            return None
        else:
            return self.db.get_round_end_time(self.round_idx)
    
    def add_to_db(self):
        """Add the round to the database, on the top of the stack of rounds."""
        self._round_idx = self.db.add_round(self)
    
    def commit(self):
        """Commit the round to the database. The round becomes immutable after this, unless reset."""
        if not self.status == "labeled":
            raise RoundStatusError(f"Cannot commit round with status {self.status}.")
        self.status = "complete"

    def reset(self):
        """Reset the round to its initial state.
        
        This will reset the status to "ready", and delete the putative library and library for experiment.
        """
        if self.status == "unknown":
            raise RoundStatusError("Cannot reset round that has not been committed to database.")
        self.db.reset_round(self.round_idx)
        self.putative_library = None
        self.library_for_exp = None
    
    def _setup_library_generation(self) -> LibraryGenerator:
        raise NotImplementedError

    def library_generation(self):
        """Generate a library of variants, store in database.
        
        Round is ready for selection after this.
        """
        # run library generator to get library, and add to db
        # set round produced for those variants to this round
        # if they were not already set
        if self.status == "ready":
            library = self.library_generator.generate()
            if len(library) == 0:
                raise EmptyLibraryException("Library generator returned empty library.")
            # this new library is connected to the database so when we change 
            # attributes of the variants, it will update the database
            library = self.db.save_library(library) 

            for variant in library.variants:
                # this may already be set, but move it up to the most recent round
                variant.round_putative = self.round_idx  # for database variants, this will update the database
                # only if this were generated from scratch and we have never seen it before do
                # we also add
                if variant.round_added is None:
                    variant.round_added = self.round_idx
            self.putative_library = library
            self.status = "generated"
            self.start_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        elif self.status == "generated":
            raise RoundStatusWarning(f"Library already generated for round {self.round_idx}, skipping. If you would like to restart this round, call round.reset()")
        else:
            raise RoundStatusError(f"Cannot generate library for round with status {self.status}")

    def _setup_library_selection(self) -> AquisitionFunction:
        raise NotImplementedError

    def library_selection(self):
        """Select a library from the putative library, and return the library for experiment."""
        if self.status == "generated":
            if self.putative_library is None:
                self.putative_library = self._attempt_recover_putative_library()
            library_for_exp = self.aquisition_function.select_library(self.putative_library)
        elif self.status == "selected":
            library_for_exp = self._attempt_recover_library_for_exp()
            raise RoundStatusWarning(f"Library already selected for round {self.round_idx}, retrieving from database. If you would like to restart this round, call round.reset()")
        else:
            raise RoundStatusError(f"Cannot select library for round with status {self.status}")
        self.library_for_exp = library_for_exp
        for variant in library_for_exp.variants:
            if variant.round_labeled is not None or :
                raise RoundStatusError(f"Cannot select library for round {self.round_idx} because variant {variant.id} was already selected for round {variant.round_selected}.")
                variant.round_selected = self.round_idx
            else:
                pass

        self.status = "selected"
        return library_for_exp
    
    def set_labels(self, labels: dict, force: bool=False, notes: str=None):
        """Set labels for the library for experiment.
        
        Params
        ------
        labels: dict
            Dictionary mapping variant id to labels.
        force: bool
            If True, overwrite existing labels.
        notes: str
            Any notes for the round.
        """
        if self.status == "selected":
            pass
        elif self.status == "labeled" and not force:
            raise RoundStatusWarning(f"Labels already set for round {self.round_idx}, skipping. If you would like to restart this round, call round.reset(). If you would like to relabel, set force=True.")
        elif self.status == "labeled" and force:
            pass
        else:
            raise RoundStatusError(f"Cannot set labels for round with status {self.status}")
        
        for id, label in labels.items():
            variant = self.library_for_exp[id]
            variant.label = label
            variant.round_labeled = self.round_idx
        self.status = "labeled"
        self.end_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.notes = notes
    
    def _attempt_recover_putative_library(self):
        """Attempt to recover putative library from database.

        This is useful if the putative library was created and committed but selection was not run.
        """
        library = self.db.get_library(round_putative=self.round_idx)
        if len(library) == 0:
            raise EmptyLibraryException(f"Attempted to recover putative library for round {self.round_idx} from database, but library was empty.")
        return library
    
    def _attempt_recover_library_for_exp(self):
        """Attempt to recover library for experiment from database.

        This is useful if the library for experiment was created and committed but lost before the library was labeled.
        """
        library = self.db.get_library(round_experiment=self.round_idx)
        if len(library) == 0:
            raise EmptyLibraryException(f"Attempted to recover library for experiment for round {self.round_idx} from database, but library was empty.")
        return library

    @classmethod
    def from_db(cls, database: CampaignDatabase, round_idx: int) -> BaseRound:
        """Load a round from the database, given the database and round index.
        
        Params
        ------
        database: CampaignDatabase
            Database to store and retrieve campaign data.
        round_idx: int
            Index of the round in the database.
        """
        params = database.get_round_params(round_idx)
        params = cls.param_class(**params)

        round = cls(database=database, params=params)
        round._round_idx = round_idx
        return round
    
"""To modify the generation or selection strategy for a round,
create a new mixin class that overloads the _setup_library_generation and
mixin to a new class that inherits from BaseRound. NoAideRound is an example
of a round that uses the default generation and selection strategies, eg. no AI
help, effectively does nothing to the library.
"""

class GenerationRoundMixin:
    def _setup_library_generation(self) -> LibraryGenerator:
        """Set library generator to take all unlabeled variants from the database.
        
        Overload to init a different library generator using round params.
        """
        from aide.zoo.library_generators import AllUnlabeledGenerator
        return  AllUnlabeledGenerator(
            database=self.db,
            round_idx=None,
        )
    
class SelectionRoundMixin:
    def _setup_library_selection(self) -> AquisitionFunction:
        """Set library selection to take all variants incoming putative library.
        
        Overload to init a different selection strategy using round params.
        """
        from aide.zoo.aquisition_functions import AllAquisition
        self.aquisition_function = AllAquisition()

class NoAideRound(GenerationRoundMixin, SelectionRoundMixin, BaseRound):
    """Minimal viable round of the campaign, where all unlabeled variants from the database are selected for experiment."""
    param_class = RoundParams