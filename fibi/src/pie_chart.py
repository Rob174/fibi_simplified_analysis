"""
Contains all the logic to build the pie chart:
- get the cases
"""
import textwrap
import pathlib as p
from string import Template
from typing import *
from dataclasses import dataclass

class CaseInputDict(TypedDict):
    avg: Literal["avgConclOk", "avgConclKo"]
    pvalue_category: Literal["small", "medium", "big"]
    es: Literal["small", "big", "nan"]

subcase_letter = Literal["a1","a2","b1","b2","c1","c2"]
case_letter = Literal["a","b","c"]

def get_case(
    average_sign_follow_hansen: bool, 
    pvalue_category: Literal["small", "medium", "big"], 
    effect_size_category: Literal["small", "big", "nan"]
) -> subcase_letter:
    case_letter = ""
    if average_sign_follow_hansen and (pvalue_category == 'small') and effect_size_category != 'small':
        case_letter = "a1"
    elif average_sign_follow_hansen and (pvalue_category == 'small') and effect_size_category == 'small':
        case_letter = "a2"
    elif not average_sign_follow_hansen and (pvalue_category == 'small') and effect_size_category != 'small':
        case_letter = "b1"
    elif not average_sign_follow_hansen and pvalue_category == 'small' and effect_size_category == 'small':
        case_letter = "b2"
    elif pvalue_category == 'big':
        case_letter = "c1"
    elif pvalue_category == 'nan':
        case_letter = "c2"
    else:
        raise Exception(f"Case not found with {average_sign_follow_hansen=} {pvalue_category=} {effect_size_category=}")
    return case_letter

@dataclass
class SubCaseACounts:
    """To store the counts of A subcases"""
    a1: int = 0
    a2: int = 0
    def __getitem__(self, item):
        return getattr(self, item)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    @property
    def total(self):
        return self.a1 + self.a2
    
@dataclass
class SubCaseBCounts:
    """To store the counts of B subcases"""
    b1: int = 0
    b2: int = 0
    b3: int = 0
    def __getitem__(self, item):
        return getattr(self, item)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    @property
    def total(self):
        return self.b1 + self.b2 + self.b3
    
@dataclass
class SubCaseCCounts:
    """To store the counts of C subcases"""
    c1: int = 0
    c2: int = 0
    def __getitem__(self, item):
        return getattr(self, item)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    @property
    def total(self):
        return self.c1 + self.c2
    
def count_cases(list_cases: List[subcase_letter]) -> Tuple[SubCaseACounts, SubCaseBCounts, SubCaseCCounts]:
    """Given a list of cases for each instances, count the number of cases in each category"""
    subcases_counters = (SubCaseACounts(), SubCaseBCounts(), SubCaseCCounts())
    for case in list_cases:
        found = False
        for subcase_counter in subcases_counters:
            if case in vars(subcase_counter):
                subcase_counter[case] += 1
                found = True
                break
        if not found:
            raise Exception
        
    return subcases_counters

class CountDict(TypedDict):
    """Contains for each category the number of samples in this category and the subcategory"""
    a: Tuple[int, Dict[str,int]]
    b: Tuple[int, Dict[str,int]]
    c: Tuple[int, Dict[str,int]]

def make_reference_str(problem: str, dataset: str, init: str = ""):
    """Builds the reference for the latex (sub)figure"""
    return f"{problem}{dataset}{init}"

def merge_and_put_tabulation(texts: Union[List[str],str], n_tabs: int = 1) -> str:
    """Merge lines of texts provided as list or string adding n_tabs tabulation"""
    if isinstance(texts, str):
        texts = texts.split("\n")
    l_texts = []
    for t in texts:
        l_texts.extend(t.split("\n"))
    tabs = "    "*n_tabs
    return ("\n" + tabs).join(l_texts)

def make_one_latex_piechart(
    problem: str, 
    dataset: str, 
    init: str, 
    subcases_counts: Tuple[SubCaseACounts, SubCaseBCounts, SubCaseCCounts], 
    template_folder: p.Path,
    main_category_text_show: bool = False,
):
    """Make a latex sunburst pie chart from the values provided:
    # In
        problem: str, to put in the title(s)
        dataset: str, to put in the title(s)
        init: str, to put in the title(s)
        subcases_counts: Tuple[SubCaseACounts, SubCaseBCounts, SubCaseCCounts], counts of each category
        template_folder: p.Path, folder where to read the templates of each line
        main_category_text_show: bool = False, wether to show the text of the main categories (A,B and C)
    # Out
        the latex code to show the sunburst pie chart
    """
    counts: CountDict = {}
    for main_case, subcase_count in zip(case_letter.__args__, subcases_counts):
        counts[main_case] = (
            sum(vars(subcase_count).values()),
            vars(subcase_count)
        )
    tot = sum(case_count for (case_count,_) in counts.values())
    def generator_inner(data: CountDict):
        angle = 0
        for main_category,(nbVals,dico_inner) in data.items():
            frac = nbVals/tot
            perc = frac*100
            start_angle = angle
            end_angle = start_angle + frac*360
            middle = (start_angle+end_angle)/2
            angle = end_angle
            yield main_category, nbVals, frac,perc, start_angle, end_angle, middle,dico_inner
    
    def generator_outer(data: CountDict):
        angle = 0
        for main_category, _, frac,perc, start_angle, end_angle, middle,dico_inner in generator_inner(data):
            for sub_category, nbValsSubCat in dico_inner.items():
                frac = nbValsSubCat/tot
                perc = frac*100
                start_angle = angle
                end_angle = start_angle + frac*360
                middle = (start_angle+end_angle)/2
                angle = end_angle
                yield main_category,sub_category, nbValsSubCat, frac,perc, start_angle, end_angle, middle
    
    with open(template_folder / "filled_arc.tex") as fp:
        template_filled_arc = Template(fp.read())
    # Outer filled arc
    background_subcategories = []
    for main_category, sub_category, _, _, perc, start_angle, end_angle, middle in generator_outer(counts):
        background_subcategories.append(
            template_filled_arc.substitute(
                start_angle=start_angle,
                end_angle=end_angle,
                pie_type="Outer",
                category=main_category
            )
        )
    # Inner filled arc
    background_categories = []
    for inner_category, _, _, perc, start_angle, end_angle, middle,dico_inner in generator_inner(counts):
        background_categories.append(
            template_filled_arc.substitute(
                start_angle=start_angle,
                end_angle=end_angle,
                pie_type="Inner",
                category=inner_category
            )
        )
    
    with open(template_folder / "text_category.tex") as fp:
        template_text_category = Template(fp.read())
    with open(template_folder / "text_percentage.tex") as fp:
        template_text_percentage = Template(fp.read())
    # Outer text
    text_subcategories = []
    for main_category, sub_category, _, _, perc, start_angle, end_angle, middle_angle in generator_outer(counts):
        text_subcategories.append(
            template_text_category.substitute(
                pie_type="Outer",
                angle=middle_angle,
                text=f"{sub_category.upper()}",
            )
        )
        if main_category_text_show or perc < 0.01:
            text_subcategories[-1] = "% "+text_subcategories[-1]
        text_subcategories.append(
            template_text_percentage.substitute(
                pie_type="Outer",
                angle=middle_angle,
                text_percentage=f"{perc:.2f}\\%",
            )
        )
        if perc < 0.01:
            text_subcategories[-1] = "% "+text_subcategories[-1]
        
    text_categories = []
    for main_category, _, _, perc, start_angle, end_angle, middle, dico in generator_inner(counts):
        text_categories.append(
            template_text_category.substitute(
                pie_type="Inner",
                angle=middle_angle,
                text=f"{main_category.upper()}",
            )
        )
        text_categories.append(
            template_text_percentage.substitute(
                pie_type="Inner",
                angle=middle_angle,
                text_percentage=f"{perc:.2f}\\%",
            )
        )
        if perc < 0.01:
            text_categories[-1] = "% "+text_categories[-1]
            text_categories[-2] = "% "+text_categories[-2]
    
        
    with open(template_folder / "one_sunburst.tex") as f:
        template_one_sunburst = Template(f.read())
    with open(template_folder / "title_subfigures.tex") as f:
        template_title_subfigures = Template(f.read())
        
    title = template_title_subfigures.substitute(
        problem=problem,
        dataset=dataset,
        init=init,
    )
    reefrence_str = make_reference_str(problem, dataset, init)
    result = template_one_sunburst.substitute(
        background_subcategories=merge_and_put_tabulation(background_subcategories, n_tabs=2),
        background_categories=merge_and_put_tabulation(background_categories, n_tabs=2),
        text_subcategories=merge_and_put_tabulation(text_subcategories, n_tabs=2),
        text_categories=merge_and_put_tabulation(text_categories, n_tabs=2),
        title=title,
        reference=reefrence_str
    )
    return result

def make_latex_piecharts_figure_for_datasets_of_problem(
    list_cases_per_init: Dict[str,List[subcase_letter]], 
    dataset: str, 
    problem: str, 
    template_folder: Union[p.Path, str]
) -> str:
    """Main entry point: to make the full figure fo a set of initializations:
    
    # Arguments
        - list_cases_per_init: Dict[str,List[subcase_letter]], for each initialization name (key) the list of subcases codes
        - dataset: str, dataset name for the title
        - problem: str, problem name for the title
        - template_folder: Union[p.Path, str], where to find the latex templates
    # Return
        - a string of the latex code to insert in an overleaf file
    """
    diagrams = {}
    for initialization, Lcases in list_cases_per_init.items():
        subcases_counts = count_cases(Lcases)
        total = sum(s.total for s in subcases_counts)
        assert total == len(Lcases), f"Total number of cases {total=} should be the same as {len(Lcases)=}"
        diagrams[initialization] = make_one_latex_piechart(
            problem=problem,
            dataset=dataset,
            init=initialization,
            subcases_counts=subcases_counts,
            template_folder=template_folder
        )
    def sort_fn(s: str) -> Tuple:
        pos = 0,0,0
        if 'RAND' in s:
            return pos
        return 1, len(s), s
    diagrams = [diagrams[k] for k in sorted(diagrams,key=sort_fn)]
    with open(template_folder / "title_figure.tex") as f:
        template_title = Template(f.read())
    with open(template_folder / "sunbursts_figure.tex") as f:
        template_figure = Template(f.read())
    result = template_figure.substitute(
        diagram=merge_and_put_tabulation(diagrams),
        title=template_title.substitute(problem=problem,dataset=dataset),
        reference=make_reference_str(problem, dataset),
    )
    return result