

from plotly.graph_objects import Figure as _GoFigure


from ..io.parse_config import config as _config
from ..core.bio_references import _BioReference
from ..core.bio_regions import BioRegion as _BioRegion
from ..core.bio_references import GenomicReference
from ..core.bio_references import TranscriptReference
from ..core.bio_references import ProteinReference
from ..core.bio_regions import BioRegion
from ..core.conversion import OutOfBoundsException as _OutOfBoundsException
from ..core.regions import Region
from ..io.gtf import GtfCluster
from ..core.regions import combine_regions as _combine_regions
from ..core.regions import get_gap as _get_gap
from ..io.gtf import collapse_genes as _collapse_genes
#from .plot_bio_regions import PlotBioRegions

import warnings as _warnings
from collections import Counter as _Counter
import plotly.graph_objects as _go
import plotly.io as _pio
import plotly.graph_objs as _pgo
Shape = _go.layout.Shape


class MissingStyleAttributeException(Exception): pass

class Figure(object):

    def __init__(self, reference: _BioReference, layout=None, style_dict=None):
        
        self.style_dict = style_dict
        
        if layout is None:
            layout = dict()
        layout.setdefault('template', 'simple_white')
        layout.setdefault('legend', dict(traceorder='reversed'))
        #layout.setdefault('hovermode', 'x') 
        layout.setdefault('margin', dict(t=100,
                                         b=80,
                                         l=30,
                                         r=80,
                                         pad=0)
                         )
        layout.setdefault('yaxis', dict(visible=False,
                                        showticklabels=False)
                         )
        
        if isinstance(reference, GenomicReference):
            xaxis_title = f'Nucelotide on chromosome {reference.chromosome}'
        elif isinstance(reference, TranscriptReference):
            xaxis_title = f'Nucelotide in transcribed region'
        elif isinstance(reference, ProteinReference):
            xaxis_title = f'Amino acid in protein'
        else:
            raise TypeError('Reference has unsupported type.')
        layout.setdefault('xaxis_title', xaxis_title)
                                  
        self.reference = reference
        self.figure = layout # call setter with layout
        self._next_update = dict(layout=layout,
                                 traces=[],
                                 shapes=[],
                                 annotations=[]
                                 )
                          
        self._x_min = None
        self._x_max = None
        self._y_min = 0.0
        self._y_max = 0.0
        self._left_buffer = 0.0
        self._right_buffer = 0.0
        

    @property
    def style_dict(self):
        return self._style_dict

    @style_dict.setter
    def style_dict(self, value):
        default_style_dict = _config['visualization']['figure_style_dict']
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

    def export(self, path, width=1200, heightPerY=100):
        """Export figure as html or common figure files with consistent look."""
        height = (self._y_max * heightPerY 
                  + self.figure.layout.margin['t']
                  + self.figure.layout.margin['b']
                  )
        if path.split('.')[-1] == 'html':
            if width == None: width = f'{width}px' #'100%'
            _pio.write_html(file=path,
                            fig=self.figure,
                            default_width=width,
                            default_height=f'{height}px')
        else:
            if width == None: width = width
            self.figure.write_image(path,
                                    width=width,
                                    height=height)
                          
    def _add_to_next_update(self, shapes=[], traces=[], layout={}, annotation=None) -> None:
        self._next_update['shapes'].extend(shapes)
        self._next_update['traces'].extend(traces)
        if annotation is not None:
            self._next_update['annotations'].append(annotation)
        for key in layout:
            self._next_update['layout'][key] = layout[key]
        return
                          
    def update_figure(self) -> None:
        self._next_update['layout']['shapes'] = self._next_update['shapes']
        self._next_update['layout']['annotations'] = self._next_update['annotations']
        self.figure.update_layout(self._next_update['layout'])
        self.figure.add_traces(self._next_update['traces'])
        self._next_update['traces'] = []
        for trace in self._figure['data']: 
            if(trace['name'] == ''): trace['showlegend'] = False
                                                  
    def _retrieve_style(self, attribute, function, style_dict):
        '''
        Retrieve the appropriate styling attribute from
        the style_dict provided to the function, or if the attribute is not
        part of the provided style_dict, retrieve it from the object style_dict
        '''
        if style_dict is not None and style_dict.get(attribute) is not None:
            return style_dict[attribute]
        elif self._style_dict.get(function) is not None and self._style_dict[function].get(attribute) is not None:
            return self._style_dict[function][attribute]
        else:
            raise KeyError(f'{attribute} not in function style_dict or config style_dict.')
               
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
        self._y_min += whitespace
        self._y_max += whitespace

    def add_regions(self, regions: list, height_ratio: float, fill_color: str, line_color=None, legend_entry=False, line_width=0.5, opacity=1.0, show_labels=False, edge_extension=0.5, row='next', row_height=1.0) -> None:
        """
        Description
        ---
        Add a list of regions to the figure in either the current, or the next row
        If the given regions reference property does not match the figure.reference,
        tries to convert to the figure reference.
        Returns true if regions where added to the plot.

        Parameters
        ---
        show_labels : False / str
                      Do not show labels on hovering over region if False. 
                      If String is provided, displays string on hover.
        """
        
        if len(regions) == 0:
            return
        
        # convert given regions to figure.reference
        converted_regions = []
        for r in regions:
            if not isinstance(r.reference, type(self.reference)):
                if not r.reference.convertible(self.reference):
                    _warnings.warn(f'Cannot plot non-convertible region {r}.')
                else:
                    try:
                        converted_regions.append(r.convert(self.reference))
                    except _OutOfBoundsException:
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
           traces.append(_go.Scatter(mode='markers',
                                      x=[None],
                                      y=[None],
                                      name=legend_entry,
                                      marker=dict(symbol='square',
                                                  opacity=opacity,
                                                  color=fill_color,
                                                  line=dict(color=line_color,
                                                            width=line_width)
                                                 )
                                      )
                         )

        # add shape for each region
        for r in regions:
            shapes.append(Shape(type='rect',
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
                try:
                    label = r.label
                except AttributeError:
                    label = ''
                if label != '':
                    hover_labels.append(f'<b>{show_labels[i]} in {r.start}-{r.end}</b><br>{label}')
                elif label == '':
                    hover_labels.append(f'<b>{show_labels[i]} in {r.start}-{r.end}</b>')                   
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
        traces = []
        if line_color is None:
            line_color = fill_color
        traces.append(_go.Scatter(mode='markers',
                                    x=[None],
                                    y=[None],
                                    name=text,
                                    marker=dict(symbol='square',
                                                opacity=opacity,
                                                color=fill_color,
                                                line=dict(color=line_color,
                                                        width=line_width)
                                                )
                                    )
                        )
        self._add_to_next_update(traces=traces)

    def add_annotation(self, text: str, side: str, margin=0) -> None:
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

        annotation = _go.layout.Annotation(text=text,
                                           #font={'size': 5},
                                           xref='paper',
                                           xanchor=x_anchor,
                                           x=x,
                                           yref='y',
                                           yanchor='middle',
                                           y=self._y_min+((self._y_max-self._y_min)/2),
                                           showarrow=False,
                                           bgcolor='rgba(255,255,255,0.95)'
                                           )
        self._add_to_next_update(annotation=annotation)

    def add_gtf_transcript_features(self, gtf_cluster, style_dict=None, transcript_ids=[], gene_ids=[], regulatory_sequence=True) -> None:
        """
        Description
        ---
        Plot the features within the given gtf_cluster. 
        If gene_id is given, only plots transcripts with the matching GtfRegion.gene_id value.
        If transcript_id is given, only plots transcripts with the matching GtfRegion.transcript_id value.
        If regulatory_sequence is True, plots regulatory sequences up- and downstream of each transcript.
        """
        
        gtf_regions = gtf_cluster.all_regions

        if len(gene_ids) != 0:
            gtf_regions = [r for r in gtf_regions if r.gene_id in gene_ids]
        if len(transcript_ids) != 0:
            gtf_regions = [r for r in gtf_regions if r.transcript_id in transcript_ids]

        transcript_ids = set([r.transcript_id for r in gtf_regions])
        try:
            transcript_ids.remove('')
        except KeyError:
            pass
                                        
        color = self._retrieve_style('color', 'add_gtf_transcript_features', style_dict)
        row_height = self._retrieve_style('row_height', 'add_gtf_transcript_features', style_dict)
        
        # Pick apart transcript regions
        transcripts = dict()
        for t in transcript_ids:
            transcript_regions = [r for r in gtf_regions if r.transcript_id==t]
            features = dict(transcript=[],
                            five_prime_utr=[],
                            start_codon=[],
                            CDS=[],
                            stop_codon=[],
                            three_prime_utr=[],
                            exon=[]
                           )
            for r in transcript_regions:
                if not isinstance(r, BioRegion):
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

        for features in minus_t+plus_t:
            plotted = []
            plotted.append(self.add_regions(regions=features['transcript'],
                                height_ratio=self._retrieve_style(f'transcript_height_ratio', 'add_gtf_transcript_features', style_dict),
                                fill_color='rgba(0,0,0,1.0)',
                                row='next',
                                row_height=row_height,
                                legend_entry=False,
                                show_labels=False,
                                line_color=color
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
                                    line_color=color
                                    ))
            for f_type in ['CDS',
                           'three_prime_utr','five_prime_utr',
                           'start_codon','stop_codon',]:                          
                plotted.append(self.add_regions(regions=features[f_type],
                                height_ratio=self._retrieve_style(f'{f_type}_height_ratio', 'add_gtf_transcript_features', style_dict),
                                fill_color=color,
                                row='current',
                                row_height=row_height,
                                legend_entry=False,
                                show_labels=f_type
                                ))
            if len(plotted) > 0 and any(plotted):
                annotation = f'{features["transcript"][0].transcript_id}, {features["transcript"][0].transcript_biotype}'
                self.add_annotation(text=annotation,
                                    side='left',
                                    margin=220)

    def add_gtf_collapsed_genes(self, gtf_cluster, style_dict=None, transcript_ids=[], gene_ids=[], regulatory_sequence=True, strand='both') -> None:
        """
        Description
        ---

        ADD OPTION TO VISUALIZE ALL STARTCODONS
        If strand is either + or -, only plots genes that are located on that strand.
        """

        gtf_regions = gtf_cluster.all_regions

        if len(gene_ids) != 0:
            gtf_regions = [r for r in gtf_regions if r.gene_id in gene_ids]
        if len(transcript_ids) != 0:
            gtf_regions = [r for r in gtf_regions if r.transcript_id in transcript_ids]
        
        if strand != 'both':
            if strand not in ['+','-']:
                raise ValueError('Invalid strand value.')
            gtf_regions = [r for r in gtf_regions if r.reference.strand == strand]

        collapsed_genes = _collapse_genes(gtf_regions)
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
                annotation = f'{g}, ({strand}), collapsed'
                self.add_annotation(text=annotation,
                                    side='left',
                                    margin=220)

            
    def add_mutations(self, regions: list, base=None, colors=None, style_dict=None, legend_entry=True, opacity=1.0, show_labels=True, row='next'):
        """
        Description
        
        Given a list of regions, selects groups of regions with the same 
        region.variant_type attribute. Each region is split into regions with length 1
        and then displayed as stacked bars, colored by its variant_type.
        """
        height_ratio = self._retrieve_style('height_ratio', 'add_mutations', style_dict)
        row_height = self._retrieve_style('row_height', 'add_mutations', style_dict)
        whitespace = row_height*(1-height_ratio)
        row_height = row_height*height_ratio
        self.add_whitespace(whitespace/2)
        row_height = self._init_row(row=row, row_height=row_height)
        
        all_mapped_locations = []
        variant_types = dict()
        for r in regions:
            locations = r.split_to_locations()
            mapped_locations = []
            unmapped_locations = []
            if isinstance(r.reference, type(self.reference)):
                mapped_locations.extend(locations)
            else:
                for l in locations:
                    try:
                        mapped_locations.append(l.convert(self.reference))
                    except _OutOfBoundsException:
                        unmapped_locations.append(l.start)
            all_mapped_locations.extend([m.start for m in mapped_locations])
            if not variant_types.get(r.variant_type):
                variant_types[r.variant_type] = []
            variant_types[r.variant_type].extend(locations)

        if len(unmapped_locations) > 0:
            _warnings.warn(f'{len(unmapped_locations)} locations could not be converted to the Figure.reference.')

        max_n_loc, max_n = _Counter(all_mapped_locations).most_common(1)[0]

        self._adjust_x(x_min=min(all_mapped_locations),
                       x_max=max(all_mapped_locations))

        shapes = []
        traces = []

        # Line at base of stacked bars
        if base is None:
            base = Region(start=min(all_mapped_locations), 
                          end=max(all_mapped_locations))
        shapes.append(Shape(type='rect',
                            x0=base.start,
                            x1=base.end,
                            y0=self._y_min,
                            y1=self._y_min,
                            line_color='rgba(0,0,0,1.0)',
                            line_width=0.5,
                            fillcolor='rgba(0,0,0,1.0)'
                            )
                      )

        stacks_height = {}
        for variant_type in sorted(variant_types.keys()):

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
            for location in variant_types[variant_type]:
                if stacks_height.get(location.start):
                    label_y_pos[location.start] = self._y_min+(stacks_height[location.start]/max_n) * (self._y_max-self._y_min) * height_ratio
                else: 
                    label_y_pos[location.start] = self._y_min

            for location in variant_types[variant_type]:

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

    def add_gtf_cluster(self, gtf_cluster: GtfCluster, strands=['-', '+'], style_dict=None, gene_ids=[], transcript_ids=[]):

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
                                        side='left',
                                        margin=220)
                    row='current'    

    def add_pas(self, pas: list, style_dict=None, legend_entry=False, annotation=True, row='next') -> None:
        """Function to plot polyadenylation and cleavage sites from the io.pas_atlas module."""

        cleavage_sites = set()
        connections = []
        for p in pas:
            cleavage_sites.add(p.cleavage_site)
            if not p.touches(p.cleavage_site):
                this_connection = _get_gap(p, p.cleavage_site)
                connections.append(_BioRegion(this_connection.start,
                                              this_connection.end,
                                              p.reference))
        connections = _combine_regions(connections)

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
        if annotation is True:
            self.add_annotation(text='Polyadenylation & cleavage signals',
                                side='left',
                                margin=220)

    def add_miRNA_binding(self, miRNA_binding_regions: list, style_dict=None, legend_entry=False, annotation=True, row='next'):
        """Function to plot miRNA binding regions from the io.targetscan module."""

        self.add_regions(regions=miRNA_binding_regions,
                         fill_color=self._retrieve_style('color', 'add_miRNA_binding', style_dict),
                         row_height=self._retrieve_style('row_height', 'add_miRNA_binding', style_dict),
                         height_ratio=self._retrieve_style('height_ratio', 'add_miRNA_binding', style_dict),
                         show_labels='miRNA-binding',
                         row=row,
                         legend_entry=legend_entry)
        if annotation is True:
            self.add_annotation(text='miRNA-binding regions',
                                side='left',
                                margin=220)

    def add_rbp_binding(self, rbp_binding_regions: list, style_dict=None, legend_entry=False, annotation=True, row='next'):

        self.add_regions(regions=rbp_binding_regions,
                         fill_color=self._retrieve_style('color', 'add_rbp_binding', style_dict),
                         row_height=self._retrieve_style('row_height', 'add_rbp_binding', style_dict),
                         height_ratio=self._retrieve_style('height_ratio', 'add_rbp_binding', style_dict),
                         show_labels='RBP-binding',
                         row=row,
                         legend_entry=legend_entry)
        if annotation is True:
            self.add_annotation(text='RBP-binding regions',
                                side='left',
                                margin=220)
        
