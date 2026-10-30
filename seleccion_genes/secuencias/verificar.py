"""Verify FASTA sources and conservative design intervals; Python >=3.10, standard library only."""
from pathlib import Path, PurePosixPath
import argparse
import collections
import csv
import hashlib
import itertools
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

EXPECTED = ('DM_ETH', 'DM_BURS_A', 'DM_ATP5F1B', 'DM_RPN6', 'DM_VATP_B',
            'EPILAB_BicC_Dm', 'EPILAB_ATPaseD_Dm', 'EPILAB_MP_MRFV', 'EPILAB_CP_MRFV')
CODONS = dict(zip((''.join(x) for x in itertools.product('TCAG', repeat=3)),
                 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))

def require(ok, message):
    if not ok:
        raise ValueError(message)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def fasta(path, strict=True):
    data = path.read_bytes()
    require(not data.startswith(b'\xef\xbb\xbf'), f'{path.name}: BOM no permitido')
    text = data.decode('ascii')
    if strict:
        require(b'\r' not in data and data.endswith(b'\n'), f'{path.name}: se esperan saltos LF y salto final')
    records, rid, lines = {}, None, []
    def finish():
        if rid is not None:
            sequence = ''.join(lines)
            require(sequence and rid not in records, f'{path.name}: registro vacío o ID duplicado {rid}')
            records[rid] = sequence
    for line in text.splitlines():
        if line.startswith('>'):
            finish()
            rid = line[1:] if strict else line[1:].split()[0]
            require(re.fullmatch(r'[A-Za-z0-9_.-]+', rid), f'{path.name}: ID inválido {rid!r}')
            lines = []
        else:
            require(rid is not None and bool(line), f'{path.name}: línea vacía o secuencia antes del encabezado')
            if strict:
                require(len(line) <= 80 and re.fullmatch('[ACGT]+', line), f'{path.name}: bases/formato inválidos en {rid}')
            lines.append(line)
    finish()
    require(records, f'{path.name}: FASTA vacío')
    return records

def embl(path):
    text = path.read_text(encoding='ascii')
    accession = re.search(r'^AC\s+(\w+);', text, re.M)[1]
    version = re.search(r'^ID\s+[^;]+; SV (\d+);', text, re.M)[1]
    sequence = re.sub(r'[0-9\s]', '', text.split('\nSQ ', 1)[1].split('\n', 1)[1].split('//', 1)[0]).upper()
    require(re.fullmatch('[ACGT]+', sequence), f'{path.name}: secuencia EMBL inválida')
    declared_md5 = re.search(r'^DR\s+MD5; ([0-9a-f]{32})\.', text, re.M)
    require(declared_md5 is not None and hashlib.md5(sequence.encode()).hexdigest() == declared_md5[1], f'{path.name}: MD5 de la secuencia fuente')
    translation = re.search(r'/translation="(.*?)"', text, re.S)[1]
    translation = re.sub(r'^FT\s+', '', translation, flags=re.M)
    translation = re.sub(r'\s', '', translation)
    coords = re.search(r'^FT\s+CDS\s+(\d+)\.\.(\d+)', text, re.M)
    return accession+'.'+version, sequence, translation, (int(coords[1]), int(coords[2]))

def read_xlsx(path):
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with zipfile.ZipFile(path) as archive:
        shared = []
        if 'xl/sharedStrings.xml' in archive.namelist():
            root = ET.fromstring(archive.read('xl/sharedStrings.xml'))
            shared = [''.join(t.text or '' for t in item.findall('.//m:t', ns)) for item in root]
        root = ET.fromstring(archive.read('xl/worksheets/sheet1.xml'))
        cells = {}
        for cell in root.findall('.//m:sheetData/m:row/m:c', ns):
            require(cell.find('m:f', ns) is None, f'Planilla fuente: fórmula inesperada en {cell.attrib["r"]}')
            if cell.attrib.get('t') == 'inlineStr':
                value = ''.join(t.text or '' for t in cell.findall('.//m:t', ns))
            else:
                element = cell.find('m:v', ns)
                value = element.text if element is not None else ''
                if cell.attrib.get('t') == 's':
                    value = shared[int(value)]
            cells[cell.attrib['r']] = value
    return cells

def gtf_cds(path, transcript, sequence):
    features = []
    for line in path.read_text(encoding='ascii').splitlines():
        if line.startswith('#'):
            continue
        cols = line.split('\t')
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', cols[8]))
        if attrs.get('transcript_id') == transcript.split('.')[0]:
            require(attrs['transcript_version'] == transcript.split('.')[1], 'Versión del transcripto GTF')
            features.append((cols, attrs))
    exons = sorted((x for x in features if x[0][2] == 'exon'), key=lambda x: int(x[1]['exon_number']))
    offset, covered = 0, []
    for exon, attrs in exons:
        lo, hi = int(exon[3]), int(exon[4])
        for feature, _ in features:
            a, b = int(feature[3]), int(feature[4])
            if feature[2] in ('CDS', 'stop_codon') and lo <= a <= b <= hi:
                if exon[6] == '+':
                    covered.extend(range(offset + a-lo+1, offset + b-lo+2))
                else:
                    covered.extend(range(offset + hi-b+1, offset + hi-a+2))
        offset += hi-lo+1
    require(offset == len(sequence), f'{transcript}: longitud/exones')
    require(covered and len(covered) == len(set(covered)), f'{transcript}: CDS vacía o solapada')
    covered.sort()
    require(covered == list(range(covered[0], covered[-1]+1)), f'{transcript}: CDS no contigua en cDNA')
    return covered[0], covered[-1]


def verify(root):
    root=root.resolve()
    sums={}
    for line in (root/'SHA256SUMS.txt').read_text(encoding='ascii').splitlines():
        expected,rel=line.split('  ',1);pure=PurePosixPath(rel)
        require(not pure.is_absolute() and '..' not in pure.parts and ':' not in rel and '\\' not in rel,'Ruta insegura')
        path=root.joinpath(*pure.parts)
        require(path.is_file() and not path.is_symlink(),f'Archivo faltante o enlace: {rel}')
        require(rel not in sums and sha(path.read_bytes())==expected,f'Integridad incorrecta: {rel}')
        sums[rel]=expected
    actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file()}
    require(actual==set(sums)|{'SHA256SUMS.txt'},'Hay archivos adicionales o faltantes')
    fixed={'README.md','genes_9.fasta','metadatos_genes.tsv','verificar.py',
        'concordancia/tramos_concordantes.fasta','concordancia/mapa_segmentos.tsv','concordancia/criterios.tsv',
        'comparacion_EPILAB/fragmentos_originales.fasta','comparacion_EPILAB/coordenadas.tsv',
        'auditoria/comparaciones.json','auditoria/diferencias.tsv','auditoria/saltos_de_alineamiento.tsv',
        'fuentes/PZ371507.1.embl','fuentes/PZ371508.1.embl','fuentes/AF265566.1.embl',
        'fuentes/Ensembl_modelos_originales.txt','fuentes/Ensembl_proteinas_originales.txt','fuentes/Ensembl_anotaciones.gtf',
        'fuentes/Andrada_transcriptos_suplemento.txt','fuentes/Andrada_Figure_S1.pdf',
        'fuentes/Sequences_EPILAB_original.xlsx','fuentes/procedencia.json'}
    require(set(sums)==fixed,'Conjunto de archivos distinto del definido en esta entrega')
    full=fasta(root/'genes_9.fasta')
    require(tuple(full)==EXPECTED,'Nueve fuentes/orden incorrectos')
    require(sum(map(len,full.values()))==12784 and len(set(full.values()))==9,'Longitud o duplicacion de las nueve secuencias')
    with (root/'metadatos_genes.tsv').open(encoding='utf-8',newline='') as f:
        primary_rows=list(csv.DictReader(f,delimiter='\t'))
    require(tuple(x['id_gen'] for x in primary_rows)==EXPECTED,'Metadatos de la entrada principal')
    for row in primary_rows:
        seq=full[row['id_gen']]
        require(len(seq)>=200 and len(seq)==int(row['longitud_nt']) and sha(seq.encode())==row['sha256_secuencia'],'Longitud/hash de entrada principal')
        rc=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
        require(all(other==row['id_gen'] or otherseq!=rc for other,otherseq in full.items()),'Duplicado por complemento inverso en entrada principal')
    segments=fasta(root/'concordancia/tramos_concordantes.fasta')
    require(len(segments)==12 and len(set(segments.values()))==12,'Se esperaban 12 tramos distintos')
    with (root/'concordancia/criterios.tsv').open(encoding='utf-8',newline='') as f:
        genes=list(csv.DictReader(f,delimiter='\t'))
    with (root/'concordancia/mapa_segmentos.tsv').open(encoding='utf-8',newline='') as f:
        mapping=list(csv.DictReader(f,delimiter='\t'))
    require(tuple(r['id_gen'] for r in genes)==EXPECTED,'Metadatos de genes incorrectos')
    require(tuple(segments)==tuple(r['id_segmento'] for r in mapping),'Orden/IDs de tramos y mapa distintos')
    models=fasta(root/'fuentes/Ensembl_modelos_originales.txt',strict=False)
    peptides=fasta(root/'fuentes/Ensembl_proteinas_originales.txt',strict=False)
    supplements=fasta(root/'fuentes/Andrada_transcriptos_suplemento.txt',strict=False)
    targets=dict(models)
    targets['Andrada_Burs']=supplements['D.maidis_burs']
    targets['Andrada_ETH_RC']=supplements['D.maidis_eth'].translate(str.maketrans('ACGT','TGCA'))[::-1]
    acc,vseq,_,_=embl(root/'fuentes/AF265566.1.embl')
    require(acc=='AF265566.1','Version viral incorrecta');targets[acc]=vseq
    source_file_text=(root/'fuentes/AF265566.1.embl').read_text(encoding='ascii')
    viral_cds=[(int(a)-1,int(b)) for a,b in re.findall(r'^FT\s+CDS\s+(\d+)\.\.(\d+)',source_file_text,re.M)]
    proof=json.loads((root/'auditoria/comparaciones.json').read_text(encoding='utf-8'))
    require(proof['package_id']=='INVESTIGATHON_FASTA_2026-10-27_DATOS','Identificador del conjunto incorrecto')
    require(tuple(proof['genes'])==EXPECTED,'Conjunto de pruebas incompleto')
    protein_ids=proof['ensembl_protein_ids']
    for tid,sequence in models.items():
        a,b=gtf_cds(root/'fuentes/Ensembl_anotaciones.gtf',tid,sequence)
        coding=sequence[a-1:b]
        require(len(coding)%3==0,f'{tid}: longitud CDS')
        peptide=''.join(CODONS[coding[i:i+3]] for i in range(0,len(coding),3))
        require(peptide.startswith('M') and peptide.endswith('*') and '*' not in peptide[:-1],f'{tid}: traduccion')
        require(peptide[:-1]==peptides[protein_ids[tid]].removesuffix('*'),f'{tid}: proteina distinta de fuente')
    cells=read_xlsx(root/'fuentes/Sequences_EPILAB_original.xlsx')
    epilab_rows={2:'EPILAB_MP_MRFV',3:'EPILAB_CP_MRFV',4:'EPILAB_BicC_Dm',5:'EPILAB_ATPaseD_Dm'}
    for number,rid in epilab_rows.items():
        raw=re.sub(r'[ \t\r\n\v\f]','',cells[f'D{number}']).upper()
        require('EPILAB_'+cells[f'B{number}']==rid and full[rid]==raw,f'{rid}: fuente EPILAB alterada')
    for row in genes:
        rid=row['id_gen'];original=full[rid];p=proof['genes'][rid]
        main_row=next(x for x in primary_rows if x['id_gen']==rid)
        require(main_row['accesion_fuente']==row['accesion_fuente'],'Accesion de entrada principal y auditoria diferentes')
        require(len(original)==int(row['longitud_fuente_nt']) and sha(original.encode())==row['sha256_fuente'],f'{rid}: longitud/hash fuente')
        a,b=int(row['consulta_inicio_fuente_0based']),int(row['consulta_fin_fuente_0based_exclusivo'])
        require(0<=a<b<=len(original),'Consulta fuera de la fuente')
        query=original[a:b]
        require(len(query)==int(row['longitud_consulta_nt']),f'{rid}: longitud consulta')
        accession=row['accesion_fuente']
        if accession:
            require(int(main_row['CDS_inicio_1based'])==a+1 and int(main_row['CDS_fin_1based_inclusivo'])==b,'Coordenadas CDS de entrada principal')
        else:
            require(main_row['CDS_inicio_1based']=='' and main_row['CDS_fin_1based_inclusivo']=='','No atribuir CDS completa no establecida a EPILAB')
        if accession.startswith('PZ'):
            public_id,public,peptide,cds=embl(root/'fuentes'/(accession+'.embl'))
            require(public_id==accession and public==original and cds==(a+1,b),f'{rid}: deposito/CDS')
            require(''.join(CODONS[query[i:i+3]] for i in range(0,len(query),3))==peptide+'*',f'{rid}: traduccion deposito')
        elif accession.startswith('ENS'):
            require(models[accession]==original and gtf_cds(root/'fuentes/Ensembl_anotaciones.gtf',accession,original)==(a+1,b),f'{rid}: modelo/CDS')
        allowed=[True]*len(query);links=[True]*(len(query)-1)
        require(p['comparisons'],f'{rid}: sin cotejo registrado')
        for c in p['comparisons']:
            target=targets[c['target_id']];lo,hi=c['target_coding_0based'];positions=c['mapping_0based']
            require(len(positions)==len(query) and len(c['base_ok'])==len(query),f'{rid}: mapa de cotejo')
            if c['target_id'].startswith('ENS'):
                g0,g1=gtf_cds(root/'fuentes/Ensembl_anotaciones.gtf',c['target_id'],target)
                require((lo,hi)==(g0-1,g1),'CDS de comparador no coincide con GTF')
            elif c['target_id']=='AF265566.1':
                require((lo,hi)==viral_cds[1 if rid=='EPILAB_MP_MRFV' else 0],'CDS viral no coincide con fuente')
            else:
                require((lo,hi)==(0,len(target)),'Contexto de suplemento inesperado')
            for q,t in enumerate(positions):
                if t is not None:require(0<=t<len(target),'Coordenada de comparador invalida')
                independently_ok=t is not None and lo<=t<hi and query[q]==target[t]
                require(bool(c['base_ok'][q])==independently_ok,'Mascara no coincide con las bases o coordenadas')
                allowed[q]&=independently_ok
            for q in range(len(links)):
                # A valid adjacent pair has unique map positions, so no insertion can be hidden.
                links[q]&=positions[q] is not None and positions[q+1] is not None and positions[q+1]==positions[q]+1
        expected_blocks=[];start=None
        for q,valid in enumerate(allowed):
            if start is not None and (not valid or (q and not links[q-1])):
                if q-start>=200:expected_blocks.append((start,q))
                start=None
            if valid and start is None:start=q
        if start is not None and len(query)-start>=200:expected_blocks.append((start,len(query)))
        listed=[r for r in mapping if r['id_gen']==rid]
        actual_blocks=[(int(r['inicio_fuente_0based'])-a,int(r['fin_fuente_0based_exclusivo'])-a) for r in listed]
        require(actual_blocks==expected_blocks,f'{rid}: tramos que incluyen posiciones no habilitadas o falta un tramo')
        require(int(row['tramos_habilitados'])==len(listed),f'{rid}: numero de tramos')
        count=0
        for m,(start,end) in zip(listed,expected_blocks):
            sid=m['id_segmento'];sequence=segments[sid];lo_source,hi_source=a+start,a+end
            require(int(m['inicio_fuente_1based'])==lo_source+1 and int(m['fin_fuente_1based_inclusivo'])==hi_source,'Conversion de coordenadas')
            require(sequence==original[lo_source:hi_source] and len(sequence)==int(m['longitud_nt'])>=200,f'{sid}: extraccion incorrecta')
            require(sha(sequence.encode())==m['sha256_secuencia'],f'{sid}: hash de bases')
            n=len(sequence)-199
            require(int(m['ventanas_200'])==n,f'{sid}: cuenta de ventanas');count+=n
            for c in p['comparisons']:
                target=targets[c['target_id']];lo,hi=c['target_coding_0based']
                # Independent string check: each entire segment is contiguous and identical.
                require(sequence in target[lo:hi],f'{sid}: falta coincidencia exacta dentro de CDS/contexto comparado')
                for i in range(n):
                    require(sequence[i:i+200] in target[lo:hi],f'{sid}: ventana de 200 discrepante')
        require(count==int(row['ventanas_200_habilitadas']),f'{rid}: total de ventanas')
    require(sum(int(r['ventanas_200']) for r in mapping)==6153,'Cuenta total inesperada')
    require(sum(map(len,segments.values()))==8541,'Longitud total de tramos inesperada')
    for sid,seq in segments.items():
        rc=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
        require(all(sid==other or otherseq!=rc for other,otherseq in segments.items()),'Duplicado por complemento inverso')
    fragments=fasta(root/'comparacion_EPILAB/fragmentos_originales.fasta')
    with (root/'comparacion_EPILAB/coordenadas.tsv').open(encoding='utf-8',newline='') as stream:
        fragment_rows=list(csv.DictReader(stream,delimiter='\t'))
    require(len(fragment_rows)==4 and len(fragments)==4,'Se esperaban cuatro controles EPILAB')
    fragment_status=[]
    for r in fragment_rows:
        rid=r['id_gen'];number=int(r['fila_fuente']);fragment=fragments[r['id_fragmento']]
        raw=re.sub(r'[ \t\r\n\v\f]','',cells[f'F{number}']).upper()
        a,b=int(r['inicio_1based_inclusivo'])-1,int(r['fin_1based_inclusivo'])
        require(fragment==raw==full[rid][a:b],f'{rid}: control de fragmento original')
        require(int(r['inicio_0based'])==a and int(r['fin_0based_exclusivo'])==b,'Coordenadas de controles')
        require(int(float(cells[f'G{number}']))==int(r['inicio_declarado_planilla']) and int(float(cells[f'H{number}']))==int(r['fin_declarado_planilla']),'Valores declarados cambiados')
        allowed_starts=set()
        for s in mapping:
            if s['id_gen']==rid:
                lo,hi=int(s['inicio_fuente_0based']),int(s['fin_fuente_0based_exclusivo'])
                allowed_starts.update(range(lo,hi-199))
        within=[x for x in range(a,b-199) if x in allowed_starts]
        require(len(within)==0,'Cambio en el estado de los controles originales que requiere revisar el informe')
        fragment_status.append({'id':rid,'ventanas_200_que_cumplen_el_filtro_conservador':len(within)})
    return {'estado':'PASS','genes':9,'nt_entrada_principal':12784,'ventanas_200_entrada_principal_por_hebra':10993,
        'tramos_opcionales':12,'nt_en_tramos_opcionales':8541,'ventanas_200_concordantes_por_hebra':6153,
        'bases_no_ACGT':0,'segmentos_identicos_a_fuentes_y_comparadores':12,
        'modelos_Ensembl_con_CDS_y_proteina_verificadas':len(models),'depositos_CDS_verificados':2,
        'secuencias_EPILAB_originales_intactas':4,'controles_EPILAB':fragment_status,'archivos_con_hash':len(sums)}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description='Verifica formato, fuentes, traducciones, tramos y concordancia de las secuencias de genes. Sin dependencias externas.')
    parser.add_argument('--json',action='store_true')
    args=parser.parse_args()
    try:result=verify(Path(__file__).resolve().parent)
    except (ValueError,KeyError,IndexError,OSError,UnicodeError,ET.ParseError,zipfile.BadZipFile) as exc:
        print(f'FAIL: {exc}',file=sys.stderr);sys.exit(1)
    if args.json:print(json.dumps(result,ensure_ascii=False,indent=2))
    else:
        print('PASS: entrada principal con 9 secuencias y 12784 nt, identicas a las fuentes declaradas.')
        print('Analisis opcional: 12 tramos concordantes, 8541 nt y 6153 ventanas de 200 nt por hebra.')
        print('Fuentes originales, CDS/proteinas y controles EPILAB comprobados.')
        print('Se verifica el conjunto de referencia; no se certifica eficacia biologica ni todas las variantes de una poblacion.')
