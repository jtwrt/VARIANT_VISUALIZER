

from plotly.graph_objects import Figure as _GoFigure


from .._config import config
from .. import core, clusters, protein_annotation
from .. import genome_annotation as ga

import warnings as _warnings
from collections import Counter as _Counter
import plotly.graph_objects as _go
import plotly.io as _pio
import plotly.graph_objs as _pgo
Shape = _go.layout.Shape
from copy import deepcopy


class MissingStyleAttributeException(Exception): 
    """Raised if required value in style_dict/config for plotting is missing"""
    pass

class Figure(object):

    def __init__(
            self, reference_type:str, cluster:clusters.Cluster, 
            transcript_ids:list='all', gene_ids:list='all',
            layout=None, style_dict=None,
            reference=None):
        """
        Description
        ---
        Class that is used to construct plots of genomic, transcript or 
        protein BioRegions.

        When initializing a Figure, provide information on which genes and transcripts
        will be visualized, as well as the reference-type of the plot.

        Visualizing multiple transcripts or genes is only possible in a Figure with
        'genomic' reference. The method add_gtf_transcript_features() is invalid for non-genomic
        figures if multiple transcripts were selected.
        The method add_gtf_collapsed_transcript_features() is only ever valid in
        genomic figures. 

        Implements methods for visualization of variants (add_variants()) 
        and different cis-regulatory elements. (add_pas(), add_miRNA_binding(), 
        add_rbp_binding()).
        
        Parameters
        ---
        layout : dict = None
            Used for initialization of the plotly.graph_objects.Figure if provided.
            Values in this dict overwrite the default values used in this class.
        style_dict : dict = None
            Overwrites config.visualization values for the instance.
            Only deined values given are overwritten.
        reference : variant_visualizer.core._BioReference
            Either GenomicReference, TranscriptReference or ProteinReference.
            All BioRegions to be displayed in this plot must be mappable
            to the chosen reference. References can best be retrieved by loading the
            correct variant_visualizer.cluster.Cluster instance and using its 
            .get_reference() method.
            The correct cluster can be identified by loading the index 
            (variant_visualizer.cluster.load_index()) and using its query methods.
        """

        self.reference_type = reference_type
        self.style_dict = style_dict
        if layout is None:
            layout = dict()
        layout.setdefault('template', 'simple_white')
        layout.setdefault('legend', dict(traceorder='reversed'))
        #layout.setdefault('hovermode', 'x') 
        layout.setdefault(
            'margin', 
            dict(
                t=20,
                b=80,
                l=50,
                r=50,
                pad=0
            )
        )
        layout.setdefault(
            'yaxis', 
            dict(
                visible=False,
                showticklabels=False
            )
         )

        self.redefine_cluster_values(
            cluster=cluster,
            gene_ids=gene_ids,
            transcript_ids=transcript_ids,
            _reference=reference
        )

        if isinstance(self.reference, core.GenomicReference):
            xaxis_title = f'Nucelotide on chromosome {self.reference.chromosome}'
        elif isinstance(self.reference, core.TranscriptReference):
            xaxis_title = f'Nucelotide in transcribed region'
        elif isinstance(self.reference, core.ProteinReference):
            xaxis_title = f'Amino acid in protein'
        else:
            raise TypeError('Reference has unsupported type.')
        layout.setdefault('xaxis_title', xaxis_title)
  
        self.figure = layout # call setter with layout
        self._next_update = dict(layout=layout,
                                 traces=[],
                                 shapes=[],
                                 annotations=[]
                                 )
                          
        self._x_min = None
        self._x_max = None
        self._y_min = 0.0 # y-range range of the track that is currently being changed
        self._y_max = 0.0
        self._left_buffer = 0.0 # determines how much whitespace is added to the x-axis, currently not in use
        self._right_buffer = 0.0
        

    @property
    def style_dict(self):
        return self._style_dict

    @style_dict.setter
    def style_dict(self, value):
        default_style_dict = config['visualization']['figure_style_dict']
        if value is None:
            self._style_dict = dict()
        else:
            self._style_dict = value
        
        for key in default_style_dict:
            self._style_dict.setdefault(key, default_style_dict[key])
            
    @property
    def figure(self):
        return self._figure

    @figure.setter
    def figure(self, layout):
        self._figure = _go.Figure(data=None, layout=layout, frames=None, skip_invalid=False)

    def redefine_cluster_values(self, cluster:clusters.Cluster, gene_ids:set='all', transcript_ids:set='all', _reference:core._BioReference=None):
        """
        Method to be used to plot non-convertible BioRegions to the same plot.
        The reference_type of the figure is not changed.
        
        Redefines internal Figure attributes the same way they are defined during initialization.

        Attributes that can be changed by invoking this method are:
            _cluster, cluster_id, _gtf_regions, 
            gene_ids, transcript_ids, default_strand_value, 
            reference
        """

        self._cluster = deepcopy(cluster)
        self.cluster_id = self._cluster.cluster_id
        self._gtf_regions = self._cluster._get_filter_regions(
            gene_ids=gene_ids,
            transcript_ids=transcript_ids
        )        
        self.gene_ids = set([r.gene_id for r in self._gtf_regions])
        self.transcript_ids = set([r.transcript_id for r in self._gtf_regions])
        if '' in self.transcript_ids: # transcript_id is '' for gene features
            self.transcript_ids.remove('')
        self._default_strand_value = self._get_default_strand_value()

        # Set figure reference value:
        if _reference is not None:
            _warnings.warn('Figure was initialized with \'_reference\' attribute. The given value will be used instead of retrieving a reference from the given cluster.')
            self.reference = _reference
        elif len(self.transcript_ids) + len(self.gene_ids) > 2:
            if self.reference_type == 'genomic':
                self.reference = self._cluster.get_reference(
                    reference_type=self.reference_type
                )
            else:
                raise ValueError('Choose \'genomic\' reference type if multiple genes and transcripts are to be displayed.')
        else:
            gene_id = None
            transcript_id = None
            if len(self.gene_ids) == 1:
                gene_id = list(self.gene_ids)[0]
            if len(transcript_ids) == 1:
                transcript_id = list(self.transcript_ids)[0]
            self.reference = self._cluster.get_reference(
                gene_id=gene_id,
                transcript_id=transcript_id,
                reference_type=self.reference_type
            )


    def export(self, path: str, width:int=750, height_per_y:int=100, html_width=None, include_plotlyjs=True):
        """
        Description
        ---
        Export figure as html or common figure files with consistent look.
        Filetype is determined from given path.
        
        Parameters
        ---
        path : str
        width : int
            width of exported figure in pixels
        height_per_y : int
            determines how narrow tracks will be displayed in the figure
        html_width : str = None
            can be given in percent \'50%\' and will overrule width for html exports
        include_plotlyjs : str = True
            determines how the plotly.js will be accesses for the interactive plot. 
            Value is directly given to plotly.graph_objects.Figure.write_html as include_plotlyjs value.
            Default value exports a independent html file.
        """
        height = (
            self._y_max * height_per_y 
            + self.figure.layout.margin['t']
            + self.figure.layout.margin['b']
        )
        if path.split('.')[-1] == 'html':
            if width == None: width = f'{width}px' #'100%'
            if html_width is not None: width = html_width
            self.figure.write_html(
                file=path,
                default_width=width,
                default_height=f'{height}px',
                include_plotlyjs=include_plotlyjs
                )
        else:
            if width == None: width = width
            self.figure.write_image(path,
                                    width=width,
                                    height=height)
                          
    def _add_to_next_update(self, shapes=[], traces=[], layout={}, annotation=None) -> None:
        """
        Internal function that gathers all updates to the plotly figure. 
        Speeds up the process of creating a figure by updating once per figure,
        instead of once per added track.
        """
        self._next_update['shapes'].extend(shapes)
        # Add traces if no identical trace is plotted already, prevent duplicate legend entries
        for t in traces:
            if all(t != prev_t for prev_t in self._next_update['traces']):
                self._next_update['traces'].append(t)
        if annotation is not None:
            self._next_update['annotations'].append(annotation)
        for key in layout:
            self._next_update['layout'][key] = layout[key]
        return
                          
    def update_figure(self) -> None:
        """
        Update the figure with all shapes, annotations and legend entries
        that were provided in previous commands.
        """
        self._next_update['layout']['shapes'] = self._next_update['shapes']
        self._next_update['layout']['annotations'] = self._next_update['annotations']
        self.figure.update_layout(self._next_update['layout'])
        self.figure.add_traces(self._next_update['traces'])
        self._next_update['traces'] = []
        for trace in self._figure['data']: 
            if(trace['name'] == ''): trace['showlegend'] = False
                                                  
    def _retrieve_style(self, attribute, function, style_dict):
        """
        Retrieve the appropriate styling attribute from
        the style_dict provided to the function, or if the attribute is not
        part of the provided style_dict, retrieve it from the object style_dict.
        """
        if style_dict is not None and style_dict.get(attribute) is not None:
            return style_dict[attribute]
        elif self._style_dict.get(function) is not None and self._style_dict[function].get(attribute) is not None:
            return self._style_dict[function][attribute]
        else:
            raise MissingStyleAttributeException(f'{attribute} not in function style_dict or config style_dict.')
               
    def _init_row(self, row: str, row_height: float) -> float:
                
        if row == 'next':
            self._y_min = self._y_max
            self._y_max += row_height
        elif row == 'current':
            pass
        else:
            raise ValueError(f'Invalid row selection {row}. Choose \'next\' or \'current\.')
        return self._y_max - self._y_min
            
    def _adjust_x(self, x_min: float, x_max: float, left_buffer=0.0, right_buffer=0.0):
        """Adjusts the x-axis range of the plotly Figure object."""
        self._left_buffer = max(self._left_buffer, left_buffer)
        self._right_buffer = max(self._right_buffer, right_buffer)
        
        if self._x_min == None:
            self._x_min = x_min
        else:
            self._x_min = min(x_min, self._x_min)

        if self._x_max == None:
            self._x_max = x_max
        else:
            self._x_max = max(x_max, self._x_max)
            
        dist = self._x_max - self._x_min
        self._add_to_next_update(layout={'xaxis':{'range': [self._x_min-(dist*self._left_buffer),
                                                            self._x_max+(dist*self._right_buffer)]}})
    
    def add_whitespace(self, whitespace: float) -> None:
        """Add an empty track to the figure. Provide the width as float."""
        self._y_min += whitespace
        self._y_max += whitespace

    def add_regions(
            self, regions: list, height_ratio: 0.9, fill_color='rgba(0,0,0,1.0)', 
            line_color=None, legend_entry:str=False, line_width=0.5, opacity=1.0, 
            show_labels=False, edge_extension=0.5, row='next', 
            row_height=1.0, _float_locations=False) -> bool:
        """
        Description
        ---
        Add a list of BioRegions to the figure. Returns a bool stating
        if regions where added to the Figure instance.
        This method is internally called by the more specialized methods.
        Mapping to different BioReferences is handled in this function.

        Parameters
        ---
        regions : list
            containing variant_visualizer.core.BioRegion object to be plotted
        height_ratio : float
            value between 0 and 1. determines how much space of the track
            will be taken up by the plottet regions.
        fill_color : str
            Fill color of the region shapes.
            Refer to plotly for valid input colors
        line_color : str = None
            If given, color of the region shape outline. Otherwise uses fill_color
        legend_entry : str = False
            Constructs a legend entry from the given string if given.
        line_width : float
            Width of the region shape outline
        opacity : float
            opacity of the region shape
        show_labels : bool or str
            Show hoverlabels on interactive plots. If True, tries to use region.label
            attribute for hoverlabels. If show_labels is a string, shows the value in
            each hoverlabel.
        edge_extension : float
            Determines how much bigger the displayed shape will. For example, 0.5
            displayes Region(1,2) from 0.5 to 2.5
        row : str
            \'next\' or \'current\': in which track regions should be plotted
        row_hieght : float
            track height
        """
        
        if len(regions) == 0:
            return False
        
        # convert given regions to figure.reference
        converted_regions = []
        for r in deepcopy(regions):
            if not isinstance(r.reference, type(self.reference)):
                r.reference.update(self.reference) # add transcript_region, coding_regions info to reference
                if not r.reference.convertible(self.reference):
                    _warnings.warn(f'Cannot plot non-convertible region {r}.')
                else:
                    try:
                        converted_regions.append(
                            r.convert(
                                self.reference,
                                _return_float_locations=_float_locations) # when converting to proteinreferences, uses float locations
                        )
                    except core.OutOfBoundsException:
                        pass
            else: 
                converted_regions.append(r)
        regions = converted_regions

        row_height = self._init_row(row=row, row_height=row_height)
        if len(regions) == 0:
            return False

        self._adjust_x(x_min=min([r.start for r in regions]), x_max=max([r.end for r in regions]))
        y_center = (self._y_max - self._y_min) * 0.5 + self._y_min
        
        traces=[]
        shapes=[]
                          
        region_height = row_height * height_ratio
        if line_color is None:
            line_color = fill_color
                          
        # add trace for legend entry
        if isinstance(legend_entry, str):
           traces.append(
               _go.Scatter(
                   mode='markers',
                    x=[None],
                    y=[None],
                    name=legend_entry,
                    marker=dict(
                        symbol='square',
                        opacity=opacity,
                        color=fill_color,
                        line=dict(
                            color=line_color,
                            width=line_width
                            )
                    )
                )
            )

        # add shape for each region
        for r in regions:
            shapes.append(
                Shape(
                    type='rect',
                    x0=r.start - edge_extension,
                    x1=r.end + edge_extension,
                    y0=y_center-0.5*region_height,
                    y1=y_center+0.5*region_height,
                    fillcolor=fill_color,
                    line_color=line_color,
                    line_width=line_width,
                    opacity=opacity
                )
            )
                          
        # Add hover labels
        if show_labels is not False:
            hover_labels = []
            for i, r in enumerate(regions):
                if isinstance(show_labels, list) and len(show_labels) == len(regions):
                    pass
                elif isinstance(show_labels, str):
                    show_labels = [show_labels for _ in regions]
                else:
                    raise TypeError('show_labels has unsupported type.')
                if isinstance(r.start, int):
                    start = str(r.start)
                    end = str(r.end)
                elif isinstance(r.start, float):
                    start = f'{r.start:.1f}'
                    end = f'{r.end:.1f}'
                else:
                    raise TypeError('Region start has invalid type.')
                if r.label is not None:
                    hover_labels.append(f'<b>{show_labels[i]} in {start}-{end}</b><br>{r.label}')
                else:
                    hover_labels.append(f'<b>{show_labels[i]} in {start}-{end}</b>') 
            traces.append(_go.Scatter(mode='markers',
                                        x=[r.start + 0.5 * (r.end - r.start) for r in regions],
                                        y=[y_center for _ in regions],
                                        name='',
                                        text=hover_labels,
                                        marker=dict(symbol='square',
                                                    opacity=0.0,
                                                    color=fill_color,
                                                    ),
                                        hovertemplate='%{text}'
                                        )
                            )
           
        self._add_to_next_update(shapes=shapes,
                                 traces=traces)
        return True

    def add_legend_entry(self, text: str, fill_color: str, symbol='square', line_width=0.5, opacity=1.0, line_color=None):
        
        # increase right margin
        current_margin = self._next_update['layout']['margin']['r']
        self._next_update['layout']['margin']['r'] = max(150, current_margin)
        
        traces = []
        if line_color is None:
            line_color = fill_color
        traces.append(_go.Scatter(mode='markers',
                                    x=[None],
                                    y=[None],
                                    name=text,
                                    marker=dict(symbol=symbol,
                                                opacity=opacity,
                                                color=fill_color,
                                                line=dict(color=line_color,
                                                        width=line_width)
                                                )
                                    )
                        )
        self._add_to_next_update(traces=traces)

    def add_annotation(self, text: str, side: str, margin=150, y='middle', y_anchor='middle') -> None:
        """Add a label to the plot, which will be aligned with the plot itself, not its coordinates."""
        if side == 'left':
            x_anchor = 'right'
            x = 0
            current_margin = self._next_update['layout']['margin']['l']
            self._next_update['layout']['margin']['l'] = max(margin, current_margin)
            text = text + ' '
        elif side == 'right':
            x_anchor = 'left'
            x = 1
            current_margin = self._next_update['layout']['margin']['r']
            self._next_update['layout']['margin']['r'] = max(margin, current_margin)
            text = ' ' + text
        else:
            raise ValueError(f'side can be either \'left\' or \'right\', not \'{side}\'.')

        if y == 'middle':
            y=self._y_min+((self._y_max-self._y_min)/2)
        elif y == 'bottom':
            y=self._y_min
        elif y == 'top':
            y=self._y_max
        else:
            raise ValueError('Invalid y location: choose middle/bottom/top')
            
        annotation = _go.layout.Annotation(
            text=text,
            #font={'size': 5},
            xref='paper',
            xanchor=x_anchor,
            x=x,
            yref='y',
            yanchor=y_anchor,
            y=y,
            showarrow=False,
            bgcolor='rgba(255,255,255,0.95)'
        )
        self._add_to_next_update(annotation=annotation)

    def add_gtf_transcript_features(self, style_dict=None, regulatory_sequence=False) -> None:
        """
        Description
        ---
        Plots the transcripts that were selected when the figure was initialized.
        This function should not be used to plot multiple transcripts outside of a
        genomic figure!
        
        Parameters
        ---
        style_dict : dict
            function style dict, to overwrite options from config.yml
        regulatory_sequence : bool
            Show regulatory sequence as defined in the variant_visualizer config.yml.
        """
    

        transcript_ids = set([r.transcript_id for r in self._gtf_regions])
        try:
            transcript_ids.remove('')
        except KeyError:
            pass
                                        
        color = self._retrieve_style('color', 'add_gtf_transcript_features', style_dict)
        row_height = self._retrieve_style('row_height', 'add_gtf_transcript_features', style_dict)
        
        # Pick apart transcript regions
        transcripts = dict()
        for t in transcript_ids:
            if len(transcript_ids) > 1 and self.reference_type != 'genomic':
                #_warnings.warn(
                #    'This function should not be used to plot multiple transcripts outside of a genomic figure!'
                #    )
                pass

            transcript_regions = [r for r in self._gtf_regions if r.transcript_id==t]
            features = dict(transcript=[],
                            five_prime_utr=[],
                            start_codon=[],
                            CDS=[],
                            stop_codon=[],
                            three_prime_utr=[],
                            exon=[]
                           )
            for r in transcript_regions:
                if not isinstance(r, ga.GtfFeature):
                    raise TypeError()
                if not features.get(r.feature):
                    features[r.feature] = []
                features[r.feature].append(r)
            transcripts[t] = features
        
        # Plot transcripts sorted by length and strandedness
        minus_t = sorted([transcripts[t] for t in transcript_ids if transcripts[t]['transcript'][0].reference.strand == '-'], 
                               key=lambda f: f['transcript'][0].get_length())
        plus_t = sorted([transcripts[t] for t in transcript_ids if transcripts[t]['transcript'][0].reference.strand == '+'], 
                               key=lambda f: f['transcript'][0].get_length())

        line_width = 1.0

        for features in minus_t+plus_t:
            plotted = []
            plotted.append(self.add_regions(regions=features['transcript'],
                                height_ratio=self._retrieve_style(f'transcript_height_ratio', 'add_gtf_transcript_features', style_dict),
                                fill_color=color,
                                row='next',
                                row_height=row_height,
                                legend_entry=False,
                                show_labels=False,
                                line_color=color,
                                line_width=line_width
                                ))
            if regulatory_sequence is True:
                for f_type in ['three_prime_regulatory_sequence','five_prime_regulatory_sequence']:
                    plotted.append(self.add_regions(regions=features[f_type],
                                    height_ratio=self._retrieve_style(f'{f_type}_height_ratio', 'add_gtf_transcript_features', style_dict),
                                    fill_color='rgba(0,0,0,0.0)',
                                    row='current',
                                    row_height=row_height,
                                    legend_entry=False,
                                    show_labels=f_type,
                                    line_color=color,
                                line_width=line_width
                                    ))
            for f_type in ['CDS',
                           'three_prime_utr','five_prime_utr',
                           'start_codon','stop_codon',]:
                # use non-black, semitransparent fill color for ProteinReference figures to see exon boundaries
                if isinstance(self.reference, core.ProteinReference):
                    fill_color = 'rgba(100,100,100,0.5)'
                    edge_extension = 1/3/2
                else:
                    fill_color = color
                    edge_extension = 0.5                
                plotted.append(self.add_regions(regions=features[f_type],
                                height_ratio=self._retrieve_style(f'{f_type}_height_ratio', 'add_gtf_transcript_features', style_dict),
                                fill_color=fill_color,
                                line_color=color,
                                row='current',
                                row_height=row_height,
                                legend_entry=False,
                                show_labels=f_type,
                                _float_locations=True,
                                edge_extension=edge_extension,
                                line_width=line_width
                                ))
            if len(plotted) > 0 and any(plotted):
                annotation = f'{features["transcript"][0].transcript_id},<br>({features["transcript"][0].reference.strand}), {features["transcript"][0].transcript_biotype}'
                self.add_annotation(text=annotation,
                                    side='left')

    def add_gtf_collapsed_transcript_features(self, style_dict=None, regulatory_sequence=False, strand='both',
                                uniprotkb_annotations:protein_annotation.UniprotAnnotations=None, legend_entry=True,
                                annotation='<br>collapsed', annotation_margin=150) -> None:
        """
        Description
        ---
        Plot the features of multiple transcripts of genes of the given cluster.
        Features of transcripts from identical genes are collapsed into a single track.
        The transcripts that will be plotted are based on the transcript_ids and gene_ids
        lists the Figure object was initalized with. 
        If the transcript_ids list is not empty, only collapses included transcripts.
        If the gene_ids list is not empty, only collapses included genes.

        Parameters
        ---
        cluster : cluster.Cluster
            variant_visualier.cluster.Cluster from where genomic regions are extraced
        style_dict : dict
            function style dict, to overwrite options from config.yml
        regulatory_sequence : bool
            Show regulatory sequence as defined in the variant_visualizer config.yml.
        strand : str
            \'both\',\'+\',\'-\' - select genes from which strand should be plotted
        uniprotkb_annotations : variant_visualizer.protein_annotation.UniprotAnnotations
            Used to exchange ensembl id annotations within the plot for gene names if provided
        legend_entry : bool
            if True, enters start and stop codon into the legend
        annotation : str
            String that is added to the standard \'gene_id, (strand)\' annotation
        annotation_margin : int
            Minimum width to which the left margin is increased, if annotations are added in this plot
        """

        if self.reference_type != 'genomic':
            raise NotImplementedError('This function can not be used outside of a genomic figure!')
            _warnings.warn(
                'This function should not be used outside of a genomic figure!'
                )

        gtf_regions = self._gtf_regions
        if strand != 'both':
            if strand not in ['+','-']:
                raise ValueError('Invalid strand value.')
            gtf_regions = [r for r in gtf_regions if r.reference.strand == strand]

        collapsed_genes = ga.collapse_gtf_genes(gtf_regions)
        gene_ids = collapsed_genes.keys()
              
        row_height = self._retrieve_style('row_height', 'add_gtf_collapsed_genes', style_dict)    
        
        for g in gene_ids:
            features = collapsed_genes[g]
            row='next'
            plotted=False
            if regulatory_sequence is True:
                for f_type in ['three_prime_regulatory_sequence','five_prime_regulatory_sequence']:
                    color = self._retrieve_style('color', 'add_gtf_collapsed_genes', style_dict)
                    if features.get(f_type) and len(features[f_type]) != 0:
                        self.add_regions(regions=features[f_type],
                                        height_ratio=self._retrieve_style(f'{f_type}_height_ratio', 'add_gtf_collapsed_genes', style_dict),
                                        fill_color=color,
                                        row=row,
                                        row_height=row_height,
                                        legend_entry=False,
                                        show_labels=f_type,
                                        line_color=color
                                        )
                        row = 'current'
                        plotted=True
            for f_type in ['intron','other_feature','CDS',
                           'three_prime_utr','five_prime_utr']:
                color = self._retrieve_style('color', 'add_gtf_collapsed_genes', style_dict)
                if features.get(f_type) and len(features[f_type]) != 0:
                    self.add_regions(regions=features[f_type],
                                    height_ratio=self._retrieve_style(f'{f_type}_height_ratio', 'add_gtf_collapsed_genes', style_dict),
                                    fill_color=color,
                                    row=row,
                                    row_height=row_height,
                                    legend_entry=False,
                                    show_labels=f_type,
                                    line_color=color
                                    )
                    row = 'current'
                    plotted=True
            for f_type in ['start_codon','stop_codon']:
                color = self._retrieve_style(f'{f_type}_color', 'add_gtf_collapsed_genes', style_dict)
                if features.get(f_type) and len(features[f_type]) != 0:
                    self.add_regions(regions=features[f_type],
                                    height_ratio=self._retrieve_style(f'{f_type}_height_ratio', 'add_gtf_collapsed_genes', style_dict),
                                    fill_color=color,
                                    row=row,
                                    row_height=row_height,
                                    legend_entry=False,
                                    show_labels=f_type,
                                    line_color=color,
                                    line_width=1.0
                                    )
                    row = 'current'
                    plotted=True

                    if legend_entry is True:
                        legend_text = {
                            'start_codon':'Start codon',
                            'stop_codon':'Stop codon'
                        }
                        self.add_legend_entry(
                            text=legend_text[f_type],
                            fill_color=color,
                        )

            if plotted is True:
                # Which strand is the transcript on
                break2 = False
                for f_type in features:
                    if break2: 
                        break
                    for r in features[f_type]:
                        strand = r.reference.strand
                        break2 = True
                        break
                if uniprotkb_annotations is not None:
                    gene_name = uniprotkb_annotations.get_gene_name(ensembl_gene_id=g,
                                                                    cluster=self._cluster)
                    text = f'{gene_name}, ({strand}){annotation}'
                else:
                    text = f'{g}, ({strand}){annotation}'
                self.add_annotation(text=text,
                                    side='left',
                                    margin=annotation_margin)

            
    def add_variants(self, variant_types='all', style_dict=None, legend_entry=True, annotation='', row='next', base='default', _variants=None):
        """
        Description
        ---
        Given a list of regions, selects groups of regions with the same 
        region.variant_type attribute. Each region is split into regions with length 1
        and then displayed as stacked bars, colored by its variant_type.

        Returns True if variants were mapped, False if not. Will always create a new track.

        Parameters
        ---
        variant_types : set = 'all'
            Only plots variants where variant_type is a value of the given set. Default plots all.
        style_dict : dict = None
            function style dict, to overwrite options from config.yml or the Figure style_dict
        legend_entry : bool = False
            If True, adds a legend entry for each annotation type plotted.
        annotation : str
            The track annotation.
        row : str
            Track position were regions should be plotted. Select 
            \'current\' or \'next\' track.
        base : BioRegion = 'default'
            Bottom line below displayed variants. Defaults to region in which variants are plotted.
        _variants : list
            List of variants that will be plotted. Selection from variant_types still applies.
            Overrules automatic retrieval of variants from Figure._cluster
        """
        height_ratio = self._retrieve_style('height_ratio', 'add_mutations', style_dict)
        row_height = self._retrieve_style('row_height', 'add_mutations', style_dict)
        whitespace = row_height*(1-height_ratio)
        row_height = row_height*height_ratio
        self.add_whitespace(whitespace/2)
        row_height = self._init_row(row=row, row_height=row_height)
        
        mapped_location_ints = []
        variants_dict = dict()

        # Add left plot annotation:
        self.add_annotation(text=annotation,
                            side='left')
        
        # Select variants to be plotted
        if _variants is not None:
            variants = _variants
        else:
            variants = self._cluster.get_variants(
                gene_ids=self.gene_ids,
                transcript_ids=self.transcript_ids
            )
        if variant_types != 'all':
            variants = [v for v in variants if v.variant_type in variant_types]
        variants = deepcopy(variants) # no lasting changes made to objects

        unmapped_locations = set()
        for r in variants:
            r.reference.strand = self.reference.strand # force overwrite variant strand value
            r.reference.update(self.reference) # add transcript_region, coding_regions info to reference
            locations = r.split_to_locations()
            mapped_locations = set()
            if isinstance(r.reference, type(self.reference)):
                [mapped_locations.add(l) for l in locations]
            else:
                for l in locations:
                    try:
                        mapped_locations.add(l.convert(self.reference))
                    except core.OutOfBoundsException:
                        unmapped_locations.add(l.start)
            mapped_location_ints.extend([m.start for m in mapped_locations])
            if not variants_dict.get(r.variant_type):
                variants_dict[r.variant_type] = []
            variants_dict[r.variant_type].extend(mapped_locations)

        shapes = []
        traces = []
        # Line at base of stacked bars, shows in which region variants were searched
        if base == 'default':
            if self.reference_type == 'genomic':
                base = core.BioRegion(
                    start=min([r.start for r in self._gtf_regions]),
                    end=max([r.end for r in self._gtf_regions]),
                    reference=self.reference
                )
            elif self.reference_type == 'transcript':
                base = core.BioRegion(
                    start=self.reference.transcript_region.start,
                    end=self.reference.transcript_region.end,
                    reference=self.reference.convert_reference_type('genomic')
                )
            elif self.reference_type == 'protein':
                base = core.BioRegion(
                    start=min([r.start for r in self.reference.coding_regions]),
                    end=max([r.end for r in self.reference.coding_regions]),
                    reference=self.reference.convert_reference_type('genomic')
                )
            else:
                raise ValueError('Invalid Figure reference type.')
        if not isinstance(base.reference, type(self.reference)):
            base = base.convert(self.reference)
        shapes.append(
            Shape(type='rect',
            x0=base.start,
            x1=base.end,
            y0=self._y_min,
            y1=self._y_min,
            line_color='rgba(0,0,0,1.0)',
            line_width=0.5,
            fillcolor='rgba(0,0,0,1.0)'
            )
        )

        if len(unmapped_locations) > 0:
            _warnings.warn(f'{len(unmapped_locations)} variant locations could not be converted to the Figure.reference. This is normal behavior when mapping variants located in regulatory sequences to a non-genomic figure.')        
        if len(mapped_location_ints) == 0:
            _warnings.warn('No variants were given or could be mapped.')
            self.add_whitespace(whitespace/2)
            self._add_to_next_update(
                shapes=shapes)
            return False

        max_n_loc, max_n = _Counter(mapped_location_ints).most_common(1)[0]

        self._adjust_x(x_min=min(mapped_location_ints),
                       x_max=max(mapped_location_ints))

        stacks_height = {}
        for variant_type in sorted(variants_dict.keys()):

            # Select color for entries of this type
            try:
                color = self._retrieve_style(
                    f'{variant_type}_color', 'add_mutations', style_dict)
            except KeyError:
                color = 'rgba(0,0,0,1.0)'
                _warnings.warn(f'Attribute \'{variant_type}_color\' is not part of the style_dict provided. Continuing with default value "{color}".')

            # Add dummy trace to create legend entry
            if legend_entry is True:
                traces.append(_go.Scatter(x=[None],
                                          y=[None],
                                          mode="lines",
                                          name=variant_type,
                                          marker={'color': color}
                                          ))
            
            labels = []
            label_y_pos = dict()
            for location in variants_dict[variant_type]:
                if stacks_height.get(location.start):
                    label_y_pos[location.start] = self._y_min+(stacks_height[location.start]/max_n) * (self._y_max-self._y_min) * height_ratio
                else: 
                    label_y_pos[location.start] = self._y_min

            for location in variants_dict[variant_type]:

                if stacks_height.get(location.start):
                    n_at_location = stacks_height[location.start]
                else: 
                    n_at_location = 0

                bottom = self._y_min + (n_at_location/max_n) * (self._y_max-self._y_min) * height_ratio
                shapes.append(_go.layout.Shape(type='rect',
                                               x0=location.start - 0.48,
                                               x1=location.start + 0.48,
                                               y0=bottom,
                                               y1=bottom + (1/max_n) * (self._y_max-self._y_min) * height_ratio,
                                               fillcolor=color,
                                               line_color=color,
                                               line_width=0.5,
                                               opacity=1.0
                                               )
                                )
                stacks_height[location.start] = n_at_location +1
                labels.append(f'{location.start}§{location.label}')
            
            counter = _Counter(labels)
            labels = dict()
            for key in counter:
                start = key.split('§')[0]
                label = key.split('§')[1]
                if not labels.get(start):
                    labels[start] = f'<b>{variant_type} at {start}</b><br>'
                labels[start] += f'{counter[key]} x {label}<br>'
            locations = labels.keys()
            traces.append(_go.Scatter(mode='markers',
                                    x=[str(l) for l in locations],
                                    y=[label_y_pos[int(l)] for l in locations],
                                    name='',
                                    text=[labels[l] for l in locations],
                                    marker={'symbol': 'line-ew',
                                            'color': color
                                            },
                                    hovertemplate='%{text}'
                            ))
        
        # Add annotation at location with most mutations
        if max_n_loc - self._x_min <= self._x_max - max_n_loc:
            annotation_text = f' {max_n}'
            annotation_xanchor = 'left'
            annotation_x = max_n_loc + 0.48
        else:
            annotation_text = f' {max_n}'
            annotation_xanchor = 'right'
            annotation_x = max_n_loc - 0.48

        annotation = _go.layout.Annotation(text=annotation_text,
                            #font={'size': 5},
                            #xref='paper',
                            xanchor=annotation_xanchor,
                            x=annotation_x,
                            yref='y',
                            yanchor='top',
                            y=self._y_min + (self._y_max-self._y_min) * height_ratio,
                            showarrow=False,
                            bgcolor='rgba(255,255,255,0.5)'
                            )
        self.add_whitespace(whitespace/2)
        self._add_to_next_update(shapes=shapes,
                                traces=traces,
                                annotation=annotation)
        return True



    def add_pas(self, strand='default', style_dict=None, legend_entry=False, annotation='PAS', row='next', _pas=None) -> None:
        """
        Description
        ---
        Plots Polyadenylation and cleavage signals and the annotated cleavage site cluster that overlapp 
        with Figure transcripts and genes.

        Parameters
        ---
        strand : str
            \'default\', \'both\', \'+\' or \'-\'. Only regions matching the strand value 
            will be plotted. Default plots all RBP-binding sites matching the strand value
            of any Figure transcript/gene.
        style_dict : dict = None
            function style dict, to overwrite options from config.yml or the Figure style_dict
        legend_entry : bool = False
            If True, adds a legend entry for each annotation type plotted.
        annotation : str = PAS
            Annotation of the track
        row : str
            Track position were regions should be plotted. Select 
            \'current\' or \'next\' track.
        _pas : list
            List of PAS that will be plotted. Filtering by strand still applies.
            Overrules automatic retrieval of variants from Figure._cluster
        """
        if _pas is None:
            pas = self._cluster.get_pas(
                gene_ids=self.gene_ids,
                transcript_ids=self.transcript_ids
            )
        else: 
            _warnings.warn(
                'Using the given PAS, not retrieving from cluster.'
            )
            pas = _pas

        if strand == 'default':
            strand = self._default_strand_value
        pas = self._select_strand_regions(
            strand=strand, 
            regions=pas
            )
        if len(pas) == 0:
            _warnings.warn('No PAS to be plotted.')
            return False

        cleavage_sites = set()
        connections = []
        for p in pas:
            cleavage_sites.add(p.cleavage_site)
            if not p.touches(p.cleavage_site):
                this_connection = core.get_gap(p, p.cleavage_site)
                connections.append(core.BioRegion(this_connection.start,
                                              this_connection.end,
                                              p.reference))
        connections = core.combine_regions(connections)

        row_height = self._retrieve_style('row_height', 'add_pas', style_dict)
        color = self._retrieve_style('color', 'add_pas', style_dict)
        self.add_regions(regions=pas,
                         fill_color=color,
                         row_height=row_height,
                         height_ratio=self._retrieve_style('signal_height_ratio', 'add_pas', style_dict),
                         show_labels='Polyadenylation signal',
                         row=row,
                         legend_entry=legend_entry)
        self.add_regions(regions=list(cleavage_sites),
                         fill_color=color,
                         row_height=row_height,
                         height_ratio=self._retrieve_style('cleavage_site_height_ratio', 'add_pas', style_dict),
                         show_labels='Cleavage site',
                         row='current')
        self.add_regions(regions=connections,
                         fill_color=color,
                         row_height=row_height,
                         height_ratio=self._retrieve_style('connection_height_ratio', 'add_pas', style_dict),
                         row='current')
    
        text = annotation
        if strand != 'both':
            text += f', ({strand})'
        self.add_annotation(text=text,
                            side='left')

    def add_miRNA_binding(self, strand='default', style_dict=None, legend_entry=False, annotation='miRNA-binding', row='next', _miRNA_binding_regions=None):
        """
        Description
        ---
        Plots miRNA-binding regions overlapping Figure transcripts and genes.

        Parameters
        ---
        strand : str
            \'default\', \'both\', \'+\' or \'-\'. Only regions matching the strand value 
            will be plotted. Default plots all RBP-binding sites matching the strand value
            of any Figure transcript/gene.
        style_dict : dict = None
            function style dict, to overwrite options from config.yml or the Figure style_dict
        legend_entry : bool = False
            If True, adds a legend entry for each annotation type plotted.
        annotation : str = 'miRNA-binding'
            Annotation of the track
        row : str
            Track position were regions should be plotted. Select 
            \'current\' or \'next\' track.
        _miRNA_binding_regions : list
            List of MiRNABinding regions that will be plotted. Filtering by strand still applies.
            Overrules automatic retrieval of variants from Figure._cluster
        """

        if _miRNA_binding_regions is None:
            miRNA_binding_regions = self._cluster.get_miRNA_binding(
                gene_ids=self.gene_ids,
                transcript_ids=self.transcript_ids
            )
        else: 
            _warnings.warn(
                'Using the given miRNA_binding_regions, not retrieving from cluster.'
            )
            miRNA_binding_regions = _miRNA_binding_regions

        if strand == 'default':
            strand = self._default_strand_value
        miRNA_binding_regions = self._select_strand_regions(
            strand=strand,
            regions=miRNA_binding_regions
        )
        if len(miRNA_binding_regions) == 0:
            _warnings.warn('No miRNA-binding regions to be plotted.')
            return False

        self.add_regions(regions=miRNA_binding_regions,
                         fill_color=self._retrieve_style('color', 'add_miRNA_binding', style_dict),
                         row_height=self._retrieve_style('row_height', 'add_miRNA_binding', style_dict),
                         height_ratio=self._retrieve_style('height_ratio', 'add_miRNA_binding', style_dict),
                         show_labels='miRNA-binding',
                         row=row,
                         legend_entry=legend_entry)
        
        text=annotation
        if strand != 'both':
            text += f', ({strand})'
        self.add_annotation(text=text,
                            side='left')

    def add_rbp_binding(self, strand='default', style_dict=None, legend_entry=False, annotation='RBP-binding', row='next', _rbp_binding_regions=None):
        """
        Description
        ---
        Plots RBP-binding regions overlapping Figure transcripts and genes.

        Parameters
        ---
        strand : str
            \'default\', \'both\', \'+\' or \'-\'. Only regions matching the strand value 
            will be plotted. Default plots all RBP-binding sites matching the strand value
            of any Figure transcript/gene.
        style_dict : dict = None
            function style dict, to overwrite options from config.yml or the Figure style_dict
        legend_entry : bool = False
            If True, adds a legend entry for each annotation type plotted.
        annotation : str = 'RBP-binding'
            Annotation of the track
        row : str
            Track position were regions should be plotted. Select 
            \'current\' or \'next\' track.
        _rbp_binding_regions : list
            List of RBPBinding regions that will be plotted. Filtering by strand still applies.
            Overrules automatic retrieval of variants from Figure._cluster
        """


        if _rbp_binding_regions is None:
            rbp_binding_regions = self._cluster.get_miRNA_binding(
                gene_ids=self.gene_ids,
                transcript_ids=self.transcript_ids
            )
        else: 
            _warnings.warn(
                'Using the given rbp_binding_regions, not retrieving from cluster.'
            )
            rbp_binding_regions = _rbp_binding_regions

        if strand == 'default':
            strand = self._default_strand_value
        rbp_binding_regions = self._select_strand_regions(
            strand=strand,
            regions=rbp_binding_regions
        )
        if len(rbp_binding_regions) == 0:
            _warnings.warn('No RBP-binding regions to be plotted.')
            return False

        self.add_regions(regions=rbp_binding_regions,
                         fill_color=self._retrieve_style('color', 'add_rbp_binding', style_dict),
                         row_height=self._retrieve_style('row_height', 'add_rbp_binding', style_dict),
                         height_ratio=self._retrieve_style('height_ratio', 'add_rbp_binding', style_dict),
                         show_labels='RBP-binding',
                         row=row,
                         legend_entry=legend_entry)
        
        text=annotation
        if strand != 'both':
            text += f', ({strand})'
        self.add_annotation(text=text,
                            side='left')
    
    def add_protein_annotations(self, protein_annotations: list, style_dict:dict=None, legend_entry=False, annotation=True, row='next'):
        """
        Description
        ---
        Add protein annotations. Annotations can be retrieved from a 
        variant_visualizer.protein_annotation.UniprotAnnotations() instance.

        Parameters
        ---
        protein_annotations : list
            list of variant_visualizer.protein_annotation.ProteinAnnotation objects
            to be plotted
        style_dict : dict = None
            function style dict, to overwrite options from config.yml or the Figure style_dict
        legend_entry : bool = False
            If True, adds a legend entry for each annotation type plotted.
        annotation : bool or str = True
            Determines the track annotation.
        row : str
            Track position were regions should be plotted. Select 
            \'current\' or \'next\' track.
        """

        function_name = 'add_uniprotkb_annotation'

        # sort by type
        a_dict = {a.annotation_type: [] for a in protein_annotations}
        [a_dict[a.annotation_type].append(a) for a in protein_annotations]

        row = row
        for a_type in a_dict.keys():
            try:
                fill_color = self._retrieve_style(
                    f'{a_type}_fill_color', function_name, style_dict
                )
            except MissingStyleAttributeException:
                fill_color = 'rgba(200,200,200,1.0)' # grey

            try:
                line_color = self._retrieve_style(
                    f'{a_type}_line_color', function_name, style_dict
                )
            except MissingStyleAttributeException:
                line_color = fill_color

            if legend_entry is True:
                this_legend_entry = a_type
            else: 
                this_legend_entry = False

            self.add_regions(
                regions=protein_annotations,
                height_ratio=self._retrieve_style(
                    'height_ratio', function_name, style_dict
                ),
                fill_color=fill_color,
                line_color=line_color,
                legend_entry=this_legend_entry,
                line_width=0.5,
                show_labels=a_type,
                row=row,
                row_height=self._retrieve_style(
                    'row_height', function_name, style_dict
                )
            )
            row = 'current'
        if annotation is not False:
            if isinstance(annotation, str):
                text = annotation
            else:
                text = 'UniprotKB annotation'
            self.add_annotation(
                text=text,
                side='left'
            )

    def _get_default_strand_value(self):
        """
        Returns both if figure gtf_regions include features on both strands,
        otherwise returns the strand value of the single strand.
        """
        strands = set()
        [strands.add(r.reference.strand) for r in self._gtf_regions]
        if len(strands) == 2:
            return 'both'
        else:
            return list(strands)[0]
    
    def _select_strand_regions(self, strand:str, regions:list):
        """Returns the selection of the regions matching thte strand values."""
        if strand == 'both':
            return regions
        elif strand in ['+', '-']:
            return [r for r in regions if r.reference.strand == strand]
        else:
            raise ValueError(f'Invalid strand value: {strand}')

    def _add_gtf_cluster(self, cluster: clusters.Cluster, strands=['-', '+'], style_dict=None):

        _warnings.warn('_add_gtf_cluster function is work in progress!')

        gtf_cluster = cluster.gtf_cluster
        plot_regions = gtf_cluster.cluster_segments

        for strand in strands:
            row = 'next'
            color = self._retrieve_style('color', 'add_gtf_cluster', style_dict)
            line_color = color
            line_width = self._retrieve_style('line_width', 'add_gtf_cluster', style_dict)
            for feature_type in plot_regions[strand]:
                if feature_type in ['five_prime_regulatory_sequence','three_prime_regulatory_sequence','five_and_three_prime_regulatory_sequence']:
                    fill_color = 'rgba(0,0,0,0.0)'
                    height_ratio = self._retrieve_style('regulatory_sequence_height_ratio', 'add_gtf_cluster', style_dict)
                elif feature_type in ['five_prime_utr','three_prime_utr','five_and_three_prime_utr']:
                    fill_color=color
                    height_ratio = self._retrieve_style('utr_height_ratio', 'add_gtf_cluster', style_dict)
                elif feature_type in ['start_codon','stop_codon']:
                    fill_color = self._retrieve_style(f'{feature_type}_color', 'add_gtf_cluster', style_dict)
                    height_ratio = self._retrieve_style(f'{feature_type}_height_ratio', 'add_gtf_cluster', style_dict)
                elif feature_type in ['CDS','intron','other_feature']:
                    fill_color=color
                    height_ratio = self._retrieve_style(f'{feature_type}_height_ratio', 'add_gtf_cluster', style_dict)
                else:
                    raise ValueError
                
                plotted = self.add_regions(plot_regions[strand][feature_type],
                                    height_ratio=height_ratio,
                                    fill_color=fill_color,
                                    line_width=line_width,
                                    line_color=line_color,
                                    show_labels=feature_type,
                                    row=row
                                    )
                if plotted is True:
                    self.add_annotation(text=f'\'{strand}\'-strand combined GTF-features',
                                        side='left')
                    row='current'    

    