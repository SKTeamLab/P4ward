from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, TextAreaField
from wtforms.validators import DataRequired

class SubmitForm(FlaskForm):
    job_name = StringField('Job name', validators=[DataRequired()])
    config_file = TextAreaField('Config', validators=[DataRequired()])
    submit = SubmitField('Run')