from flask import Flask
from flask import render_template, url_for, flash, redirect

from forms import SubmitForm

app = Flask(__name__)
app.config['SECRET_KEY'] = 'eb3d6960e536c85410af4127fcdf70fe'

@app.route("/")
@app.route("/home")
def home():
    return(render_template('home.html'))

@app.route('/about')
def about():
    return(render_template('about.html'))

@app.route('/submit', methods=['GET','POST'])
def submit():
    form = SubmitForm()
    if form.validate_on_submit():
        flash('Job Submitted', 'success')
        return(redirect(url_for('home')))
    return(render_template('submit.html', form=form))
