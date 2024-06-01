"""
Contains all the logic to build the pie chart:
- get the cases
"""
import pathlib as p
from string import Template
from typing import *
from dataclasses import dataclass
import fibi.src.statistical_test as statistical_test
import numpy as np

class CaseInputDict(TypedDict):
    avg: Literal["avgConclOk", "avgConclKo"]
    pvalue_category: Literal["small", "medium", "big"]
    es: Literal["small", "big", "nan"]

subcase_letter = Literal["FI","fi","BI","bi","NC"]
case_letter = Literal["f","b","n"]

def get_case(
    fi_better: bool, 
    pvalue_category: Literal["small", "medium", "big"], 
    effect_size_category: Literal["small", "big", "nan"]
) -> subcase_letter:
    case_letter = ""
    if fi_better and (pvalue_category == 'small') and effect_size_category != 'small':
        case_letter = "FI"
    elif fi_better and (pvalue_category == 'small') and effect_size_category == 'small':
        case_letter = "fi"
    elif not fi_better and (pvalue_category == 'small') and effect_size_category != 'small':
        case_letter = "BI"
    elif not fi_better and pvalue_category == 'small' and effect_size_category == 'small':
        case_letter = "bi"
    elif pvalue_category == 'big' or pvalue_category == "nan":
        case_letter = "NC"
    else:
        raise Exception(f"Case not found with {fi_better=} {pvalue_category=} {effect_size_category=}")
    return case_letter

@dataclass
class SubCaseFICounts:
    """To store the counts of A subcases"""
    FI: int = 0
    fi: int = 0
    def __getitem__(self, item):
        return getattr(self, item)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    @property
    def total(self):
        return self.FI + self.fi
    
@dataclass
class SubCaseBICounts:
    """To store the counts of B subcases"""
    BI: int = 0
    bi: int = 0
    def __getitem__(self, item):
        return getattr(self, item)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    @property
    def total(self):
        return self.bi + self.BI
    
@dataclass
class SubCaseNCCounts:
    """To store the counts of C subcases"""
    NC: int = 0
    def __getitem__(self, item):
        return getattr(self, item)
    def __setitem__(self, key, value):
        return setattr(self, key, value)
    @property
    def total(self):
        return self.NC
    
def count_cases(list_cases: List[subcase_letter]) -> Tuple[SubCaseFICounts, SubCaseBICounts, SubCaseNCCounts]:
    """Given a list of cases for each instances, count the number of cases in each category"""
    subcases_counters = (SubCaseFICounts(), SubCaseBICounts(), SubCaseNCCounts())
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
    FI: Tuple[int, Dict[str,int]]
    BI: Tuple[int, Dict[str,int]]
    NC: Tuple[int, Dict[str,int]]

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
    subcases_counts: Tuple[SubCaseFICounts, SubCaseBICounts, SubCaseNCCounts], 
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
            biggest_subcat_num = max(dico_inner.values())
            one_cat_only = (nbVals - biggest_subcat_num)/tot < 0.01
            yield main_category, nbVals, frac,perc, start_angle, end_angle, middle,dico_inner, one_cat_only
    
    def generator_outer(data: CountDict):
        angle = 0
        for main_category, _, frac,perc, start_angle, end_angle, middle,dico_inner, one_cat_only in generator_inner(data):
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
                category=main_category.upper()
            )
        )
    # Inner filled arc
    background_categories = []
    for inner_category, _, _, perc, start_angle, end_angle, middle,dico_inner, one_cat_only in generator_inner(counts):
        background_categories.append(
            template_filled_arc.substitute(
                start_angle=start_angle,
                end_angle=end_angle,
                pie_type="Inner",
                category=inner_category.upper()
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
                text=f"\\textit{{{sub_category}}}",
            )
        )
        if perc < 0.01:
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
    for main_category, _, _, perc, start_angle, end_angle, middle, dico, one_cat_only in generator_inner(counts):
        text_categories.append(
            template_text_category.substitute(
                pie_type="Inner",
                angle=middle,
                text=f"{perc:.2f}\\%",
            )
        )
        if perc < 0.01 or main_category_text_show or one_cat_only:
            text_categories[-1] = "% "+text_categories[-1]
            
    
        
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
    def sort_fn(s: str) -> int:
        if s == 'RAND' or s == "random":
            return 0
        elif s == "greedy" or s == "GREEDY":
            return 1
        else:
            return 2
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

def get_case_from_diff(diff: np.ndarray, maximization: bool = True, init_random: bool = True) -> subcase_letter:
    avg_diff = np.mean(diff)
    test_result = statistical_test.run_wilcoxon(diff)
    pvalue_category = statistical_test.mapping_pvalue_category(test_result['pvalue'])
    effect_size_category = statistical_test.mapping_effect_size_category(test_result["effect_size"], test="Wilcoxon")
    fi_better = statistical_test.fi_is_better(avg_diff, maximization=maximization)
    
    subcase_letter_str: subcase_letter = get_case(fi_better, pvalue_category, effect_size_category)
    return subcase_letter_str