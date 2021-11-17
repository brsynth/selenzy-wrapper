#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 12:20:49 2017

@author: Pablo Carbonell, jerrywzy
"""
import os, subprocess, glob, time, shutil
import argparse, uuid, json, csv
import logging
from logging.handlers import RotatingFileHandler
import Selenzy
from flask import Flask, flash, render_template, request, redirect, url_for, send_from_directory, jsonify
from flask_restful import Resource, Api
from flask import session
from werkzeug import secure_filename
import pandas as pd
import numpy as np

global session
app = Flask(__name__)
api = Api(app)
app.config['SECRET_KEY'] = str(uuid.uuid4())
app.config['MARVIN'] = False
app.config['KEEPDAYS'] = 10
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024
ALLOWED_EXTENSIONS = set(['txt', 'rxn', 'smi', 'smarts', 'smirks', 'csv', 'fasta', 'fas', 'fa'])


def arguments():
    parser = argparse.ArgumentParser(description='Options for the webserver')
    parser.add_argument('-uploaddir', default='uploads',
                        help='Upload folder')
    parser.add_argument('-datadir', default='data',
                        help='Data directory for required databases files')
    parser.add_argument('-logdir', default='log',
                        help='Logging folder')    
    parser.add_argument('-d', action='store_true',
                        help='Run in debug mode (no preload)')
    arg = parser.parse_args()
    return arg

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.',1)[1].lower() in ALLOWED_EXTENSIONS

def file_path(uniqueid, filename):
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    uniquename = os.path.join(uniquefolder, filename)
    return uniquename

def save_rxn(rxninfo):
    global session
    filename = secure_filename(rxninfo.filename)
    try:
        uniquename = file_path(session['uniqueid'], filename)
    except:
        init_session()
        uniquename = file_path(session['uniqueid'], filename)
    rxninfo.save(uniquename)
    outname = file_path(session['uniqueid'], session['uniqueid'])
    rxninfo = Selenzy.sanitizeRxn(uniquename, outname)
    session['rxninfo'] = rxninfo
    return rxninfo

def init_session():
    global session
    maintenance(app.config['KEEPDAYS'])
    reset_session()
    uniqueid = session['uniqueid']
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    if not os.path.exists(uniquefolder):
        os.mkdir(uniquefolder)
    session['uniquefolder'] = uniquefolder
    session['rxnifo'] = None
    session['status'] = False
    session['username'] = session['uniqueid']
    # Restart the Score for each new session
    session['SCORE'] = Selenzy.seqScore()


def reset_session():
    global session
    uniqueid = str(uuid.uuid4())
    app.logger.info( 'New session: %s' % (uniqueid,) )
    session['uniqueid'] = uniqueid

def run_session(rxntype, rxninfo, targets, direction, host, fp, noMSA):
    global session
    uniqueid = session['uniqueid']
    uniquefolder = session['uniquefolder']
    csvfile = "selenzy_results.csv"
    app.logger.info( 'Run session: %s' % (uniqueid,) )
    success, app.config['TABLES'] = Selenzy.analyse(['-'+rxntype, rxninfo], 
                                                    targets,
                                                    app.config['DATA_FOLDER'],  
                                                    uniquefolder,
                                                    csvfile,
                                                    pdir = int(direction),
                                                    host = host,
                                                    fp = fp,
                                                    NoMSA = noMSA,
                                                    pc = app.config['TABLES']
    ) # this creates CSV file in Uploads directory
    if success:
        data = Selenzy.updateScore(file_path(uniqueid, csvfile), session['SCORE'])
        return data, csvfile, uniqueid

    
def retrieve_session(csvinfo):
    global session
    uniqueid = session['uniqueid']
    uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], uniqueid)
    if not os.path.exists(uniquefolder):
        os.mkdir(uniquefolder)
    filename = secure_filename(csvinfo.filename)
    uniquename = file_path(uniqueid, filename)
    csvinfo.save(uniquename)
    data = pd.read_csv(uniquename)
    data.index = data.index + 1
    csvfile = os.path.basename(uniquename)
    data.rename_axis('Select', axis="columns")
    return data, csvfile, uniqueid


def maintenance(expDay=10):
    secs = expDay*24*60*60
    for folder in glob.glob(os.path.join(app.config['UPLOAD_FOLDER'], '*')):
        name = os.path.basename(folder)
        if name.startswith('debug'):
            continue
        modiftime = os.path.getmtime(folder)
        lapse = time.time() - modiftime
        if lapse > secs:
            # Double check that this an upload folder containing reactions
            if  len( glob.glob( os.path.join(folder, '*.rxn') ) ) > 0:
                try:
                    for x in glob.glob(os.path.join(folder, '*')):
                        os.unlink(x)
                except:
                    pass
            try:
                os.rmdir(folder)
                app.logger.info( 'Clean up: %s' % (folder,) )
            except:
                pass

        
class RestGate(Resource):
    """ REST interface, returns api info """
    def get(self):
        return {'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem'}

class RestQuery(Resource):
    """ REST interface to Selenzy, by default it does not run the MSA to be faster. 
    We init an independent session for the REST request."""
    def post(self):
        global session
        args = request.json
        init_session()
        if 'rxnid' in args and 'db' in args and 'smarts' not in args:
            """ Retrieve the SMARTS from the database id """
            db = args['db']
            rxnid = db+':'+args['rxnid']
            if rxnid in app.config['TABLES'].rxnref and app.config['TABLES'].rxnref[rxnid] in app.config['TABLES'].smir:
                mnxrid = app.config['TABLES'].rxnref[rxnid]
                smarts = app.config['TABLES'].smir[ mnxrid ][0]
                if mnxrid in app.config['TABLES'].rxndir:
                    if app.config['TABLES'].rxndir[mnxrid] == '-1':
                        smarts = app.config['TABLES'].smir[ mnxrid ][1]
                try:
                    outname = file_path(session['uniqueid'], session['uniqueid'])
                    rxninfo = Selenzy.sanitizeSmarts(smarts, outname)
                    args['smarts'] = rxninfo
                except:
                    pass            

        if 'smarts' in args:
            """ Submit SMARTS query """
            rxntype = 'smarts'
            rxninfo = args['smarts']
            if 'targets' in args:
                targets = args['targets']
            else:
                targets = '50'
            if 'direction' in args:
                direction = int(args['direction'])
            else:
                direction = 0
            if 'noMSA' in args:
                noMSA = args['noMSA']
            else:
                noMSA = True
            if 'host' in args:
                host = args['host']
            else:
                host = '83333'
            if 'fp' in args:
                fp = args['fp']
            else:
                fp = 'RDK'
            if 'score' in args:
                session['SCORE'] = Selenzy.seqScore(args['score'])
            try:
                if isinstance(rxninfo, (list, tuple) ):
                    data = []
                    for instance in rxninfo:
                        dat, csvfile, sessionid = run_session(rxntype, instance, targets, direction, host, fp, noMSA)
                        data.append(dat)
                    data = pd.DataFrame(data)
                else:
                    data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, host, fp, noMSA)
                return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': data.to_json()})
            except:
                return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': None})
        else:
            return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': None})

class RestSource(Resource):
    """ REST interface, returns api info """
    def get(self):
        orgs = {}
        for seq in app.config['ORG']:
            orgs[app.config['ORG'][seq][1]] = app.config['ORG'][seq][0]
        return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': orgs})

class RestFinger(Resource):
    """ REST interface, returns api info """
    def get(self):
        fp = Selenzy.availableFingerprints()
        return jsonify({'app': 'Selenzy', 'version': '1.0', 'author': 'Synbiochem', 'data': list(fp)})


api.add_resource(RestGate, '/REST')

api.add_resource(RestQuery, '/REST/Query')

api.add_resource(RestSource, '/REST/Source')
api.add_resource(RestFinger, '/REST/Fingerprints')


@app.errorhandler(404)
def page_not_found(e):
    return redirect(url_for('upload_form'))

@app.route('/')
def upload_form():
    if 'username' not in session:
        return redirect(url_for('login'))
    return render_template("my_form.html", username=session['username'], fingerprints=Selenzy.availableFingerprints())


@app.route('/login', methods=['GET', 'POST'])
def login():
    if app.debug == True:
      session['username'] = 'debug'  
      init_session()
      return redirect(url_for('upload_form'))
    if request.method == 'POST':
        session['username'] = request.form['username']
        init_session()
        return redirect(url_for('upload_form'))
    else:
        init_session()
        return redirect(url_for('upload_form'))
    return '''
    <form method="post">
    <p><input type=text name=username>
    <p><input type=submit value=Login>
    </form>
    '''

@app.route('/logout')
def logout():
    session.pop('username', None)
    return redirect(url_for('login'))


@app.route('/msa', methods=['POST'])
def post_msa():
    """ Post safely the MSA """
    if request.method == 'POST':
        sessionid = json.loads(request.values['sessionid'])
        msafile = os.path.join(app.config['UPLOAD_FOLDER'], sessionid, 'sequences_aln.fasta')
        treefile = os.path.join(app.config['UPLOAD_FOLDER'], sessionid, 'sequences.dnd')
        if os.path.exists(msafile) and os.path.exists(treefile):
            msa = open(msafile).readlines()
            tree = open(treefile).readlines()
            return json.dumps({'msa': ''.join(msa), 'tree': ' '.join(tree)})
    return redirect ( url_for('upload_form') )

@app.route('/msaview', methods=['GET'])
def display_msa():
    """ Display the MSA """
    if request.method == 'GET':
        if 'id' in request.values:
            sessionid = request.values['id']
            msafile = os.path.join(app.config['UPLOAD_FOLDER'], sessionid, 'sequences_aln.fasta')
            if os.path.exists(msafile):
                return render_template('viewmsa.html', sessionid=sessionid)
    return redirect ( url_for('upload_form') )

@app.route('/display', methods=['POST'])
def display_reaction(marvin=app.config['MARVIN']):
    """ Validates the query and displays the reaction """
    if request.method == 'POST':
        size = (600,400)
        if 'file' in request.files and len(request.files['file'].filename) > 0:
            fileinfo = request.files['file']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            rxninfo = save_rxn(fileinfo)
            success = True
            if len(rxninfo) == 0:
                success = False
                data = ''
            else:
                if marvin:
                    svgstream = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'],
                                                         outname = str(uuid.uuid4()), marvin=True)
                    data = svgstream.decode('utf-8')
                    if len(data) == 0:
                        success = False
                else:
                    outfile, size = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'],
                                                             outname = str(uuid.uuid4()), marvin=False)
                    if len(outfile) == 0:
                        success = False
                    data = os.path.join('/results', session['uniqueid'], 'files', os.path.basename(outfile))
                session['rxninfo'] = rxninfo
                session['rxntype'] = 'smarts'
                session['status'] = True
                success = True
            return json.dumps( {'data': data, 'status': session['status'], 'success': success, 'svg': marvin, 'size': size, 'smarts': rxninfo} )
        elif len(request.form['smarts']) > 0:
            outname = file_path(session['uniqueid'], session['uniqueid'])
            rxninfo = Selenzy.sanitizeSmarts(request.form['smarts'], outname)
            success = True
            if marvin:
                svgstream = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=True)
                data = svgstream.decode('utf-8')
                if len(data) == 0:
                    success = False
            else:
                outfile, size = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=False)
                if len(outfile) == 0:
                    success = False
                data = os.path.join('/results', session['uniqueid'], 'files', os.path.basename(outfile))
            session['rxninfo'] = rxninfo
            session['rxntype'] = 'smarts'
            session['status'] = True
            return json.dumps( {'data': data, 'status': session['status'], 'success': success, 'svg': marvin, 'size': size, 'smarts': rxninfo} )
        elif len(request.form['rxnid']) > 0:
            db = request.form['rdb']
            rxnid = db+':'+request.form['rxnid']
            if rxnid in app.config['TABLES'].rxnref and app.config['TABLES'].rxnref[rxnid] in app.config['TABLES'].smir:
                mnxrid = app.config['TABLES'].rxnref[rxnid]
                smarts = app.config['TABLES'].smir[ mnxrid ][0]
                if mnxrid in app.config['TABLES'].rxndir:
                    if app.config['TABLES'].rxndir[mnxrid] == '-1':
                        smarts = app.config['TABLES'].smir[ mnxrid ][1]
                outname = file_path(session['uniqueid'], session['uniqueid'])
                rxninfo = Selenzy.sanitizeSmarts(smarts, outname)
                success = True
                if marvin:
                    svgstream = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=True)
                    data = svgstream.decode('utf-8')
                    if len(data) == 0:
                        success = False
                else:
                    outfile, size = Selenzy.display_reaction(rxninfo, outfolder=session['uniquefolder'], outname = str(uuid.uuid4()), marvin=False)
                    if len(outfile) == 0:
                        success = False
                    data = os.path.join('/results', session['uniqueid'], 'files', os.path.basename(outfile))
                session['rxninfo'] = rxninfo
                session['rxntype'] = 'smarts'
                session['status'] = True
                return json.dumps( {'data': data, 'status': session['status'], 'success': success, 'svg': marvin, 'size': size, 'smarts':smarts} )

@app.route('/sorter', methods=['POST'])
def sort_table():
    """ Sorts table """
    if request.method == 'POST':
        jfilter = json.loads(request.values.get('filter'))
        try:
            filt = [int(x) for x in jfilter]
        except:
            return
        session = json.loads(request.values.get('session'))
        csvname = os.path.basename(json.loads(request.values.get('csv')))
        csvfile = os.path.join(app.config['UPLOAD_FOLDER'], session, csvname)
        outdir = os.path.join(app.config['UPLOAD_FOLDER'], session)
        head, rows = Selenzy.read_csv(csvfile)
        sortrows = Selenzy.sort_rows(rows, filt)
        Selenzy.updateMSA(outdir, sortrows)
        Selenzy.write_csv(csvfile, head, sortrows)
        data = pd.read_csv(csvfile)
        data.index = data.index + 1
        data.rename_axis('Select', axis="columns")
        return json.dumps( {'data': {'csv':  data.to_html(), 'filter': filt}} )

@app.route('/adder', methods=['POST'])
def add_rows():
    """ Add rows to table """
    if request.method == 'POST':
        if 'session' in request.values:
            sessionid = request.values['session']
        else:
            flash("Bad request")
            return redirect (request.url)            
        if 'fasta' in request.files and len(request.files['fasta'].filename) > 0:
            fileinfo = request.files['fasta']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            uniquefolder = os.path.join(app.config['UPLOAD_FOLDER'], sessionid)
            fastafile = sessionid+'.fasta'
            uniquename = os.path.join(uniquefolder, fastafile)
            fileinfo.save(uniquename)
            dndFile = os.path.join(uniquefolder, 'sequences.dnd')
            if os.path.exists(dndFile):
                noMSA = False
            else:
                noMSA = True
            csvfile = Selenzy.extend_sequences('sequences.fasta', fastafile, uniquefolder, noMSA)
            data = Selenzy.updateScore(file_path(uniquefolder, csvfile), session['SCORE'])

            data = pd.read_csv(csvfile)
            data.index = data.index + 1
            data.rename_axis('Select', axis="columns")
            # TO DO: update fasta file
            return json.dumps( {'data': {'csv':  data.to_html()}} )


@app.route('/remover', methods=['POST'])
def delete_rows():
    """ Sorts table """
    if request.method == 'POST':
        selrows = json.loads(request.values.get('filter'))
        session = json.loads(request.values.get('session'))
        csvname = os.path.basename(json.loads(request.values.get('csv')))
        outdir = os.path.join(app.config['UPLOAD_FOLDER'], session)
        csvfile = os.path.join(outdir, csvname)
        head, rows = Selenzy.read_csv(csvfile)
        filt = []
        for i in selrows:
            try:
                index = int(i) - 1
                filt.append(index)
            except:
                continue
        newtargets = []
        newrows = []
        for j in range(0, len(rows)):
            if j not in filt:
                newtargets.append(rows[j][head.index('Seq. ID')])
                newrows.append(rows[j])
        fastaFile = os.path.join(outdir, "sequences.fasta")
        Selenzy.write_fasta(fastaFile, newtargets, app.config['TABLES'])
        # Avoid issues with sequence ids
        fastaShortNameFile = os.path.join(outdir, "seqids.fasta")
        Selenzy.write_fasta(fastaShortNameFile, newtargets, app.config['TABLES'], short=True)
        # Recompute MSA if exists
        dndFile = os.path.join(outdir, 'sequences.dnd')
        if os.path.exists(dndFile):
            cons = Selenzy.doMSA(fastaShortNameFile, outdir)
            for i in range(0, len(newrows)):
                try:
                    newrows[i][head.index('Consv. Score')] = cons[newrows[i][head.index('Seq. ID')]]
                except:
                    pass
        Selenzy.write_csv(csvfile, head, newrows)
        data = pd.read_csv(csvfile)
        data.index = data.index + 1
        data.rename_axis('Select', axis="columns")
        return json.dumps( {'data': {'csv':  data.to_html()}} )
                          

@app.route('/scorer', methods=['POST'])
def score_table():
    """ Score table """
    if request.method == 'POST':
        score = json.loads(request.values.get('score'))
        sessid = json.loads(request.values.get('session'))
        csvname = os.path.basename(json.loads(request.values.get('csv')))
        csvfile = os.path.join(app.config['UPLOAD_FOLDER'], sessid, csvname)
        session['SCORE'] = Selenzy.seqScore(score)
        data = Selenzy.updateScore(csvfile, session['SCORE'])
        return json.dumps( {'data': {'csv':  data.to_html()}} )

@app.route('/debug', methods=['GET'])
def show_table():
    if app.debug == True:
        csvfile = os.path.join(app.config['UPLOAD_FOLDER'], 'debug', 'selenzy_results.csv')
        data = Selenzy.updateScore(csvfile, session['SCORE'])
        sessionid = 'debug'
        data.rename_axis('Select', axis="columns")
        return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid,
                               flags={'fasta': False, 'msa': False}, score=session['SCORE'])
    else:
        return redirect ( url_for('upload_form') )

@app.route('/results', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        """ The POST request should come from an already initalised session """
        if 'uniqueid' not in session:
            return redirect ( url_for('upload_form') )
        # check if post request has smarts part
        if 'csv' in request.files  and len(request.files['csv'].filename) > 0:
            fileinfo = request.files['csv']   
            if fileinfo.filename == '' or not allowed_file(fileinfo.filename):
                flash("No file selected")
                return redirect (request.url)
            data, csvfile, sessionid = retrieve_session(fileinfo)
            return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid, flags={'fasta': False, 'msa': False}, score=session['SCORE'])
        else:
            try:
                rxninfo = session['rxninfo']
                rxntype = session['rxntype']
            except:
                return redirect(url_for('login'))

        direction = 0
        noMSA = False
        targets = request.form['targets']
        host = request.form['host']
        fp = request.form['finger']
        if request.form.get('direction'):
            direction = 1
        if request.form.get('noMSA'):
            noMSA = True
        try:
            data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, host, fp, noMSA)
            return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid, flags={'fasta': True, 'msa': not noMSA}, score=session['SCORE'])
        except:
            return redirect( url_for("upload_form") )
    elif request.method == 'GET':
        if request.args.get('fasta') is not None:
            """ This is a request that is handled by ajax, do not return anything """
            sessionid = request.args.get('session')
            return ('', 204)
        else:
            """ A GET request would require an independently initialised session """
            init_session()
            smarts = request.args.get('smarts')
            if smarts is None:
                return redirect( url_for("upload_form") )
            host = request.args.get('host')
            if host is None:
                host = '83333'
            fp = request.args.get('fp')
            if fp is None:
                fp = 'RDK'
            rxntype = 'smarts'
            rxninfo = smarts
            direction = 0
            noMSA = False
            targets = 20
            session['rxninfo'] = rxninfo
            session['rxntype'] = rxntype
            try:
                data, csvfile, sessionid = run_session(rxntype, rxninfo, targets, direction, host, fp, noMSA)
                return render_template('results.html', tables=data.to_html(), csvfile=csvfile, sessionid=sessionid,
                                       flags={'fasta': True, 'msa': not noMSA}, score=session['SCORE'])
            except:
                return redirect( url_for("upload_form") )
    return redirect( url_for("upload_form") )
    
@app.route('/results/<sessionid>/files/<filename>')
def results_file(sessionid,filename):
    return send_from_directory(os.path.join(app.config['UPLOAD_FOLDER'], sessionid), filename)

# Reconfigure for gunicorn
if __name__== "__main__":  #only run server if file is called directly

    arg = arguments()
    

    app.config['UPLOAD_FOLDER'] = os.path.abspath(arg.uploaddir)
    app.config['LOG_FOLDER'] = os.path.abspath(arg.logdir)
    app.config['DATA_FOLDER'] = os.path.abspath(arg.datadir)

    if arg.d:
        app.config['DEBUG'] = True
        app.config['PRELOAD'] = True
    else:
        app.config['DEBUG'] = False
        app.config['PRELOAD'] = True        

    app.config['ORG'] = Selenzy.seqOrganism(arg.datadir, "seq_org.tsv")

    if app.config['PRELOAD']:
        app.config['TABLES'] = Selenzy.readData(arg.datadir)
    else:
        app.config['TABLES'] = None

    handler = RotatingFileHandler(os.path.join(app.config['LOG_FOLDER'], 'selenzy.log'), maxBytes=10000, backupCount=1)

    log = logging.getLogger('werkzeug')
    log.addHandler(handler)
    app.logger.addHandler(handler)

    app.run(host="0.0.0.0",port=5000, debug=app.config['DEBUG'], threaded=True)
#    app.run(port=5000, debug=True)
